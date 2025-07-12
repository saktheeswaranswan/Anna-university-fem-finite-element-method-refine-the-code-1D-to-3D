let numNodes = 4;
let numElements = 3;
let fem;
let sliderNodes, sliderElems, dofInput, showButton, connectivityButton, exportButton;
let highlightDOF = -1;
let showConnectivity = false;

function setup() {
  createCanvas(1280, 720);
  textFont('monospace');

  sliderNodes = createSlider(3, 10, numNodes, 1).position(20,20).input(resetFEM);
  sliderElems = createSlider(2, 9, numElements, 1).position(20,50).input(resetFEM);
  
  dofInput = createInput('').position(20,90).size(60);
  showButton = createButton('🔎 Highlight DOF').position(90,90)
    .mousePressed(() => highlightDOF = int(dofInput.value()));
  
  connectivityButton = createButton('📋 Toggle Connectivity').position(220,90)
    .mousePressed(() => showConnectivity = !showConnectivity);
  
  exportButton = createButton('💾 Export Results').position(380,90)
    .mousePressed(() => fem.exportResults());
  
  resetFEM();
}

function draw() {
  background(255);
  fem.drawMatrix(highlightDOF);
  fem.drawForceVector(highlightDOF);

  fill(0);
  text(`Nodes: ${sliderNodes.value()} | Elements: ${sliderElems.value()}`, 20, 140);
  text('Global Stiffness Matrix (Sparse view)', 320, 30);
  text('Force Vector', 1020, 30);

  if (showConnectivity) {
    fem.drawConnectivity();
  } else {
    text('Toggle "📋 Connectivity" to view element DOF mapping.', 20, 320);
  }
}

function resetFEM() {
  numNodes = sliderNodes.value();
  numElements = sliderElems.value();
  fem = new TrussFEM(numNodes, numElements);
  fem.assemble();
}

// =======================
class TrussFEM {
  constructor(NN, NE) {
    this.NN = NN;
    this.NE = NE;
    this.DOF = 2 * NN;
    this.K = Array(this.DOF).fill().map(() => Array(this.DOF).fill(0));
    this.F = Array(this.DOF).fill(0);
    this.E = 2e7;
    this.Area = 1.0;
    this.nodes = [];
    this.elements = [];
    this.SEcollection = [];
    this.generateMesh();
  }

  generateMesh() {
    this.nodes = [];
    this.elements = [];
    for (let i = 0; i < this.NN; i++) {
      let x = (i % 2 === 0) ? 0 : 40;
      let y = floor(i / 2) * 40;
      this.nodes.push([x, y]);
    }
    for (let e = 0; e < this.NE; e++) {
      let n1 = e % this.NN;
      let n2 = (e + 1) % this.NN;
      if (n1 !== n2) this.elements.push([n1, n2]);
    }
  }

  assemble() {
    this.SEcollection = [];
    this.K.forEach(row => row.fill(0));
    this.F.fill(0);

    for (let e = 0; e < this.elements.length; e++) {
      let [i1,i2] = this.elements[e];
      let [x1,y1] = this.nodes[i1];
      let [x2,y2] = this.nodes[i2];
      let dx = x2 - x1, dy = y2 - y1;
      let L = sqrt(dx*dx + dy*dy);
      let c = dx/L, s = dy/L;
      let EAL = this.E*this.Area/L;

      let SEloc = [
        [ c*c,  c*s, -c*c, -c*s],
        [ c*s,  s*s, -c*s, -s*s],
        [-c*c, -c*s,  c*c,  c*s],
        [-c*s, -s*s,  c*s,  s*s]
      ].map(r => r.map(v => v * EAL));

      let dofMap = [2*i1,2*i1+1,2*i2,2*i2+1];
      this.SEcollection.push({element:e, dofMap, matrix:SEloc});

      for (let i=0;i<4;i++) for (let j=0;j<4;j++){
        this.K[dofMap[i]][dofMap[j]] += SEloc[i][j];
      }
    }

    // apply fixed supports at node 0
    this.K[0][0] += 1e20;
    this.K[1][1] += 1e20;
    // apply some load at last node
    this.F[this.DOF-2] = 20000;
    this.F[this.DOF-1] = -25000;
  }

  drawMatrix(highlight) {
    let s=10, bx=300, by=50;
    for (let i=0;i<this.DOF;i++){
      for (let j=0;j<this.DOF;j++){
        if (abs(this.K[i][j])>1e-3){
          let blink = frameCount%60 < 30;
          fill(highlight>=0 && (i===highlight||j===highlight) && blink ? color(255,0,0) : 0);
          rect(bx+j*s, by+i*s, s, s);
        }
      }
    }
  }

  drawForceVector(highlight) {
    let s=10, bx=1020, by=50;
    for (let i=0;i<this.F.length;i++){
      if (abs(this.F[i])>1e-3){
        let blink = frameCount%60 < 30;
        fill(i===highlight && blink ? color(255,0,0) : 0);
        rect(bx, by+i*s, s, s);
      }
    }
  }

  drawConnectivity() {
    let bx=850, by=340;
    textAlign(LEFT);
    this.SEcollection.forEach(({element, dofMap}, idx) => {
      text(`Elem ${element}: DOFs [${dofMap.join(', ')}]`, bx, by + idx*20);
    });
  }

  exportResults() {
    let lines = [];
    lines.push("=== Sparse Global Stiffness Entries ===");
    for (let i = 0; i < this.DOF; i++){
      for (let j = 0; j < this.DOF; j++){
        if (Math.abs(this.K[i][j]) > 1e-3) {
          lines.push(`K[${i}][${j}] = ${this.K[i][j].toExponential(4)}`);
        }
      }
    }
    lines.push("\n=== Force Vector (non-zero) ===");
    this.F.forEach((f,i) => {
      if (Math.abs(f) > 1e-3) lines.push(`F[${i}] = ${f.toExponential(4)}`);
    });
    lines.push("\n=== Element DOF Connectivity ===");
    this.SEcollection.forEach(({element, dofMap}) => {
      lines.push(`Element ${element}: DOFs = [${dofMap.join(', ')}]`);
    });

    saveStrings(lines, 'truss_sparse_results.txt');
  }
}

