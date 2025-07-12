// === 2D TRUSS FEM STIFFNESS MATRIX + FORCE VECTOR VISUALIZER ===
// Written in p5.js | Matrix Assembly Format like TRUSS2D (Chandrupatla)

let numNodes = 4;
let numElements = 4;
let fem;
let sliderNodes, sliderElems, dofInput, showButton;
let highlightDOF = -1;

function setup() {
  createCanvas(1280, 720);
  textFont('monospace');

  sliderNodes = createSlider(3, 8, numNodes, 1);
  sliderNodes.position(20, 20);
  sliderNodes.input(() => resetFEM());

  sliderElems = createSlider(3, 12, numElements, 1);
  sliderElems.position(20, 50);
  sliderElems.input(() => resetFEM());

  dofInput = createInput('');
  dofInput.position(20, 90);
  dofInput.size(60);

  showButton = createButton('🔎 Highlight DOF');
  showButton.position(90, 90);
  showButton.mousePressed(() => highlightDOF = int(dofInput.value()));

  resetFEM();
}

function draw() {
  background(255);
  fem.drawMatrix(highlightDOF);
  fem.drawForceVector(highlightDOF);

  fill(0);
  text(`Nodes: ${sliderNodes.value()} | Elements: ${sliderElems.value()}`, 20, 140);
  text('Global Stiffness Matrix', 320, 30);
  text('Force Vector', 1020, 30);
}

function resetFEM() {
  numNodes = sliderNodes.value();
  numElements = sliderElems.value();
  fem = new TrussFEM(numNodes, numElements);
  fem.assemble();
}

// ==== CLASSES ====
class TrussFEM {
  constructor(NN, NE) {
    this.NN = NN; // nodes
    this.NE = NE; // elements
    this.DOF = NN * 2;
    this.K = Array(this.DOF).fill().map(() => Array(this.DOF).fill(0));
    this.F = Array(this.DOF).fill(0);
    this.E = 2e7;
    this.Area = 1.0;
    this.nodes = [];
    this.elements = [];
    this.generateMesh();
  }

  generateMesh() {
    for (let i = 0; i < this.NN; i++) {
      let x = (i % 2 == 0) ? 0 : 40;
      let y = floor(i / 2) * 30;
      this.nodes.push([x, y]);
    }

    for (let i = 0; i < this.NE; i++) {
      let n1 = i % this.NN;
      let n2 = (i + 1) % this.NN;
      if (n1 != n2) this.elements.push([n1, n2]);
    }
  }

  assemble() {
    for (let e = 0; e < this.elements.length; e++) {
      let [i1, i2] = this.elements[e];
      let [x1, y1] = this.nodes[i1];
      let [x2, y2] = this.nodes[i2];
      let dx = x2 - x1, dy = y2 - y1;
      let L = sqrt(dx*dx + dy*dy);
      let c = dx / L;
      let s = dy / L;
      let EAL = this.E * this.Area / L;

      // Element stiffness matrix SE (4x4)
      let SE = [
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s]
      ];

      // DOF mapping
      let dofMap = [2*i1, 2*i1+1, 2*i2, 2*i2+1];

      for (let i = 0; i < 4; i++) {
        let I = dofMap[i];
        for (let j = 0; j < 4; j++) {
          let J = dofMap[j];
          this.K[I][J] += EAL * SE[i][j];
        }
      }
    }

    // Example boundary condition + load (hardcoded)
    this.K[0][0] += 1e20; this.K[1][1] += 1e20; // Fix node 0 (x,y)
    this.F[this.DOF - 2] = 20000; // Last node x-direction force
    this.F[this.DOF - 1] = -25000; // Last node y-direction force
  }

  drawMatrix(highlight) {
    let s = 10;
    let baseX = 300, baseY = 50;
    for (let i = 0; i < this.DOF; i++) {
      for (let j = 0; j < this.DOF; j++) {
        if (abs(this.K[i][j]) > 1e-3) {
          let blink = ((frameCount % 60) < 30);
          if (highlight >= 0 && (i === highlight || j === highlight) && blink) fill(255, 0, 0);
          else fill(0);
          rect(baseX + j * s, baseY + i * s, s, s);
        }
      }
    }
  }

  drawForceVector(highlight) {
    let s = 10;
    let baseX = 1020, baseY = 50;
    for (let i = 0; i < this.F.length; i++) {
      if (abs(this.F[i]) > 1e-3) {
        let blink = ((frameCount % 60) < 30);
        if (i === highlight && blink) fill(255, 0, 0);
        else fill(0);
        rect(baseX, baseY + i * s, s, s);
      }
    }
  }
}
