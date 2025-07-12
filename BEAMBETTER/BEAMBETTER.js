let numNodes = 3;
let numElements = 2;
let fem;
let sliderNodes, sliderElems, dofInput, showButton, connectivityButton;
let highlightDOF = -1;
let showConnectivity = false;

function setup() {
  createCanvas(1280, 720);
  textFont('monospace');

  sliderNodes = createSlider(3, 20, numNodes, 1);
  sliderNodes.position(20, 20);
  sliderNodes.input(() => resetFEM());

  sliderElems = createSlider(2, 19, numElements, 1);
  sliderElems.position(20, 50);
  sliderElems.input(() => resetFEM());

  dofInput = createInput('');
  dofInput.position(20, 90);
  dofInput.size(60);

  showButton = createButton('🔎 Highlight DOF');
  showButton.position(90, 90);
  showButton.mousePressed(() => highlightDOF = int(dofInput.value()));

  connectivityButton = createButton('📋 Show Connectivity');
  connectivityButton.position(220, 90);
  connectivityButton.mousePressed(() => showConnectivity = !showConnectivity);

  resetFEM();
}

function draw() {
  background(255);
  fem.drawMatrix(highlightDOF);
  fem.drawForceVector(highlightDOF);

  fill(0);
  text(`Nodes: ${sliderNodes.value()} | Elements: ${sliderElems.value()}`, 20, 140);
  text('Global Stiffness Matrix (non-banded view)', 320, 30);
  text('Force Vector', 1020, 30);

  if (showConnectivity) {
    fem.drawElementStiffnessWithDOF();
  } else {
    text('Toggle "📋 Show Connectivity" to view per-element stiffness matrices and DOF mapping.', 20, 320);
  }
}

function resetFEM() {
  numNodes = sliderNodes.value();
  numElements = sliderElems.value();
  fem = new BeamFEM(numNodes, numElements);
  fem.assemble();
}

// ==== CLASSES ====

class BeamFEM {
  constructor(NN, NE) {
    this.NN = NN;
    this.NE = NE;
    this.DOF = NN * 2; // DOF: vertical + rotation per node
    this.K = Array(this.DOF).fill().map(() => Array(this.DOF).fill(0));
    this.F = Array(this.DOF).fill(0);
    this.E = 2e5;
    this.I = 4e6;
    this.nodes = [];
    this.elements = [];
    this.SEcollection = [];
    this.generateMesh();
  }

  generateMesh() {
    let spacing = 1000;
    for (let i = 0; i < this.NN; i++) {
      this.nodes.push([i * spacing, 0]);
    }
    for (let i = 0; i < this.NE; i++) {
      this.elements.push([i, i + 1]);
    }
  }

  assemble() {
    for (let e = 0; e < this.elements.length; e++) {
      let [i1, i2] = this.elements[e];
      let [x1, _] = this.nodes[i1];
      let [x2, __] = this.nodes[i2];
      let L = Math.abs(x2 - x1);
      let EIL = this.E * this.I / Math.pow(L, 3);

      let SE = [
        [12 * EIL, 6 * EIL * L, -12 * EIL, 6 * EIL * L],
        [6 * EIL * L, 4 * EIL * L * L, -6 * EIL * L, 2 * EIL * L * L],
        [-12 * EIL, -6 * EIL * L, 12 * EIL, -6 * EIL * L],
        [6 * EIL * L, 2 * EIL * L * L, -6 * EIL * L, 4 * EIL * L * L]
      ];

      let dofMap = [2 * i1, 2 * i1 + 1, 2 * i2, 2 * i2 + 1];
      this.SEcollection.push({ element: e, dofMap, matrix: SE });

      for (let i = 0; i < 4; i++) {
        let I = dofMap[i];
        for (let j = 0; j < 4; j++) {
          let J = dofMap[j];
          this.K[I][J] += SE[i][j];
        }
      }
    }

    // Boundary conditions
    this.K[0][0] += 1e20;
    this.K[1][1] += 1e20;
    this.F[this.DOF - 2] = -6000;
    this.F[this.DOF - 1] = 1e6;
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
          ellipse(baseX + j * s + s / 2, baseY + i * s + s / 2, s * 0.7, s * 0.7);
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
        ellipse(baseX + s / 2, baseY + i * s + s / 2, s * 0.7, s * 0.7);
      }
    }
  }

  drawElementStiffnessWithDOF() {
    let s = 7;
    let baseX = 20, baseY = 340;
    for (let idx = 0; idx < this.SEcollection.length; idx++) {
      let { element, matrix, dofMap } = this.SEcollection[idx];
      let xOffset = baseX + (idx % 4) * 300;
      let yOffset = baseY + floor(idx / 4) * 120;

      fill(0);
      text(`Element ${element} | DOF: [${dofMap.join(', ')}]`, xOffset, yOffset);
      for (let i = 0; i < 4; i++) {
        for (let j = 0; j < 4; j++) {
          let val = matrix[i][j].toExponential(2);
          fill(0);
          text(val, xOffset + j * 70, yOffset + 20 + i * 15);
        }
      }
    }
  }
}

