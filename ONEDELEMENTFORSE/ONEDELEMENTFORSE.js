// Optimized 1D FEM Visualization in p5.js with UI + Force Vector

let numElements = 5;
let highlightNode = -1;
let femSolver;

let slider, nodeInput, showButton, exportButton;
let forceBox, forceButton;
let connectivityText = "";
let forceText = "";

function setup() {
  createCanvas(1200, 720);
  textFont('monospace');

  slider = createSlider(2, 100, numElements, 1);
  slider.position(20, 20);
  slider.input(() => {
    numElements = slider.value();
    femSolver = new FEMSolver(numElements);
    femSolver.solve();
  });

  nodeInput = createInput('');
  nodeInput.position(20, 60);
  nodeInput.size(100);

  showButton = createButton('▶ Show Connectivity');
  showButton.position(130, 60);
  showButton.mousePressed(() => {
    highlightNode = int(nodeInput.value());
    connectivityText = femSolver.getConnectivityText(highlightNode);
  });

  forceBox = createInput('');
  forceBox.position(20, 100);
  forceBox.size(100);
  forceBox.attribute('placeholder', 'Node for Force');

  forceButton = createButton('⚡ Show Force');
  forceButton.position(130, 100);
  forceButton.mousePressed(() => {
    let id = int(forceBox.value());
    forceText = femSolver.getForceText(id);
  });

  exportButton = createButton('📁 Export Results');
  exportButton.position(20, height - 30);
  exportButton.mousePressed(() => femSolver.exportResults());

  femSolver = new FEMSolver(numElements);
  femSolver.solve();
}

function draw() {
  background(255);
  femSolver.drawStiffnessMatrix(highlightNode);
  femSolver.drawDisplacementPlot();

  fill(0);
  text("Displacement Plot", 20, height - 120);
  text("Stiffness Matrix (Red = connected to input node)", 270, femSolver.K.length * 10 + 40);

  fill(240);
  rect(20, 150, 250, 180);
  fill(0);
  textSize(12);
  text(connectivityText, 25, 165);
  text(forceText, 25, 280);
}

class Node {
  constructor(id, x) {
    this.id = id;
    this.x = x;
    this.u = 0;
  }
}

class Element {
  constructor(id, n1, n2, E, A) {
    this.id = id;
    this.n1 = n1;
    this.n2 = n2;
    this.Le = n2.x - n1.x;
    this.k = E * A / this.Le;
    this.strain = 0;
    this.stress = 0;
  }

  computeStressStrain(E) {
    this.strain = (this.n2.u - this.n1.u) / this.Le;
    this.stress = E * this.strain;
  }

  isConnectedTo(nodeId) {
    return this.n1.id === nodeId || this.n2.id === nodeId;
  }
}

class FEMSolver {
  constructor(nElem) {
    this.E = 210e9;
    this.A = 0.01;
    this.L = 1.0;
    this.nElem = nElem;
    this.nodes = [];
    this.elements = [];
    this.K = [];
    this.F = [];
    this.setup();
  }

  setup() {
    const nNodes = this.nElem + 1;
    const dx = this.L / this.nElem;
    this.nodes = [];
    for (let i = 0; i < nNodes; i++) {
      this.nodes.push(new Node(i, i * dx));
    }
    this.elements = [];
    for (let i = 0; i < this.nElem; i++) {
      this.elements.push(new Element(i, this.nodes[i], this.nodes[i + 1], this.E, this.A));
    }
    this.K = Array(nNodes).fill().map(() => Array(nNodes).fill(0));
    this.F = Array(nNodes).fill(0);
  }

  assemble() {
    const n = this.nodes.length;
    for (let e of this.elements) {
      const i = e.n1.id, j = e.n2.id;
      const k = e.k;
      this.K[i][i] += k;
      this.K[j][j] += k;
      this.K[i][j] -= k;
      this.K[j][i] -= k;
    }
    this.F[n - 1] = 1000.0;
    this.K[0][0] += 1e20;
    this.F[0] += 0;
  }

  solve() {
    this.assemble();
    const u = this.gaussianElimination(this.K, this.F);
    for (let i = 0; i < u.length; i++) this.nodes[i].u = u[i];
    for (let e of this.elements) e.computeStressStrain(this.E);
  }

  gaussianElimination(A, b) {
    const n = b.length;
    let M = A.map(row => row.slice());
    let f = b.slice();
    let x = Array(n).fill(0);
    for (let k = 0; k < n - 1; k++) {
      for (let i = k + 1; i < n; i++) {
        const factor = M[i][k] / M[k][k];
        for (let j = k; j < n; j++) M[i][j] -= factor * M[k][j];
        f[i] -= factor * f[k];
      }
    }
    for (let i = n - 1; i >= 0; i--) {
      let sum = 0;
      for (let j = i + 1; j < n; j++) sum += M[i][j] * x[j];
      x[i] = (f[i] - sum) / M[i][i];
    }
    return x;
  }

  drawDisplacementPlot() {
    const x0 = 250, y0 = height - 100;
    const scaleX = 500, scaleY = 1e5;
    stroke(0);
    noFill();
    beginShape();
    for (let n of this.nodes) {
      let x = x0 + n.x * scaleX;
      let y = y0 - n.u * scaleY;
      vertex(x, y);
    }
    endShape();
  }

  drawStiffnessMatrix(highlight) {
    const s = 10;
    push();
    translate(270, 20);
    for (let i = 0; i < this.K.length; i++) {
      for (let j = 0; j < this.K.length; j++) {
        if (Math.abs(this.K[i][j]) > 1e-6) {
          if (highlight >= 0 && (i === highlight || j === highlight)) fill(255, 0, 0);
          else fill(0);
          rect(j * s, i * s, s, s);
        }
      }
    }
    pop();
  }

  getConnectivityText(nodeId) {
    if (nodeId < 0 || nodeId >= this.nodes.length) return "❌ Invalid node.";
    let txt = "Node " + nodeId + " connected to:\n";
    for (let e of this.elements) {
      if (e.isConnectedTo(nodeId)) {
        txt += `Element ${e.id}: Node ${e.n1.id} ↔ Node ${e.n2.id}\n`;
      }
    }
    return txt;
  }

  getForceText(nodeId) {
    if (nodeId < 0 || nodeId >= this.F.length) return "❌ Invalid node for force.";
    return `Force at Node ${nodeId}: ${this.F[nodeId].toFixed(2)} N`;
  }

  exportResults() {
    let output = "==== Node Displacements (m) ====";
    for (let n of this.nodes) output += `\nNode ${n.id}: ${n.u}`;
    output += "\n\n==== Element Strain ====";
    for (let e of this.elements) output += `\nElement ${e.id}: ${e.strain}`;
    output += "\n\n==== Element Stress (Pa) ====";
    for (let e of this.elements) output += `\nElement ${e.id}: ${e.stress}`;
    output += "\n\n==== Global Stiffness Matrix ====";
    for (let row of this.K) output += "\n" + row.map(v => v.toFixed(2)).join("\t");

    output += "\n\n==== Global Force Vector ====";
    for (let i = 0; i < this.F.length; i++) output += `\nF[${i}] = ${this.F[i].toFixed(2)}`;

    saveStrings([output], "fem_results.txt");
    print("✅ Exported fem_results.txt");
  }
} 
