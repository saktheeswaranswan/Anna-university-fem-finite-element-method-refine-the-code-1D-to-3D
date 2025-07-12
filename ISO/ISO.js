// === 2D ISOPARAMETRIC QUAD FEM VISUALIZER ===
// 4‑node quadrilateral elements, 2×2 Gauss integration (Example 8.4)

let fem;
let dofInput, showButton;
let highlightDOF = -1;

function setup() {
  createCanvas(1280, 720);
  textFont('monospace');
  fem = new QuadFEM();
  fem.assemble();

  dofInput = createInput('').position(20, 20).size(60);
  showButton = createButton('🔎 Highlight DOF')
    .position(90, 20)
    .mousePressed(() => highlightDOF = int(dofInput.value()));
}

function draw() {
  background(255);
  fem.drawMesh();
  fem.drawMatrix(highlightDOF);
  fem.drawForceVector(highlightDOF);

  fill(0);
  text('Enter DOF# (0–17) and click to highlight K & F', 20, 60);
}

class QuadFEM {
  constructor() {
    // 9 nodes, 4 quads
    this.nodes = [
      [0,0],[30,0],[60,0],
      [0,15],[30,15],[60,15],
      [0,30],[30,30],[60,30]
    ];
    this.elements = [[0,3,4,1],[1,4,5,2],[3,6,7,4],[4,7,8,5]];
    this.DOF = this.nodes.length * 2;
    this.K = Array(this.DOF).fill().map(() => Array(this.DOF).fill(0));
    this.F = Array(this.DOF).fill(0);
    this.E = 7e4;
    this.nu = 0.33;
    this.t = 10;
  }

  assemble() {
    const gauss = [-1/Math.sqrt(3), 1/Math.sqrt(3)];
    const D = this.getDMatrix();

    this.elements.forEach(el => {
      let SE = Array(8).fill().map(() => Array(8).fill(0));
      gauss.forEach(xi => gauss.forEach(eta => {
        let {B, detJ} = this.getB(el, xi, eta);
        for (let i=0; i<8; i++) {
          for (let j=0; j<8; j++) {
            let val =
              B[0][i]*D[0][0]*B[0][j] +
              B[0][i]*D[0][1]*B[1][j] +
              B[1][i]*D[1][0]*B[0][j] +
              B[1][i]*D[1][1]*B[1][j] +
              B[2][i]*D[2][2]*B[2][j];
            SE[i][j] += val * detJ * this.t;
          }
        }
      }));

      let dof = [];
      el.forEach(n => dof.push(2*n, 2*n+1));
      for (let i=0; i<8; i++) for (let j=0; j<8; j++) {
        this.K[dof[i]][dof[j]] += SE[i][j];
      }
    });

    // fix left-edge nodes 0,1,3,4
    [0,1,3,4].forEach(n => {
      this.K[2*n][2*n] += 1e20;
      this.K[2*n+1][2*n+1] += 1e20;
    });

    // apply load at node 9 x-direction -> DOF index 16 (0-based)
    this.F[16] = -10000;
  }

  getDMatrix() {
    let E=this.E, v=this.nu;
    let c = E/(1 - v*v);
    return [[c, c*v, 0],[c*v, c, 0],[0,0,c*(1-v)/2]];
  }

  getB(el, ξ, η) {
    // shape function derivatives
    let dNξ = [-(1-η),(1-η),(1+η),-(1+η)].map(v=>v/4);
    let dNη = [-(1-ξ),-(1+ξ),(1+ξ),(1-ξ)].map(v=>v/4);

    // Jacobian
    let J = [[0,0],[0,0]];
    el.forEach((n,i) => {
      let [x,y] = this.nodes[n];
      J[0][0] += dNξ[i]*x;
      J[0][1] += dNξ[i]*y;
      J[1][0] += dNη[i]*x;
      J[1][1] += dNη[i]*y;
    });
    let detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    let invJ = [
      [ J[1][1]/detJ, -J[0][1]/detJ ],
      [ -J[1][0]/detJ, J[0][0]/detJ ]
    ];

    // B matrix
    let B = Array(3).fill().map(() => Array(8).fill(0));
    el.forEach((n,i) => {
      let dx = invJ[0][0]*dNξ[i] + invJ[0][1]*dNη[i];
      let dy = invJ[1][0]*dNξ[i] + invJ[1][1]*dNη[i];
      B[0][2*i]   = dx;
      B[1][2*i+1] = dy;
      B[2][2*i]   = dy;
      B[2][2*i+1] = dx;
    });

    return {B, detJ};
  }

  drawMesh() {
    stroke(0); fill(0);
    this.nodes.forEach((p,i) => {
      let X = 200 + p[0]*5;
      let Y = 600 - p[1]*5;
      ellipse(X,Y,5);
      text(i+1, X+5, Y-5);
    });
    noFill(); stroke(150);
    this.elements.forEach(el => {
      beginShape();
      el.forEach(n => {
        let [x,y] = this.nodes[n];
        vertex(200 + x*5, 600 - y*5);
      });
      endShape(CLOSE);
    });
  }

  drawMatrix(h) {
    let s=8, bx=600, by=50;
    for (let i=0;i<this.DOF;i++) for (let j=0;j<this.DOF;j++) {
      if (abs(this.K[i][j])>1e-2) {
        let blink = frameCount % 60 < 30;
        fill((h===i||h===j)&&blink ? 'red' : 'black');
        rect(bx+j*s, by+i*s, s, s);
      }
    }
  }

  drawForceVector(h) {
    let s=8, bx=600+this.DOF*8+20, by=50;
    for (let i=0;i<this.DOF;i++) {
      if (abs(this.F[i])>1e-2) {
        let blink = frameCount % 60 < 30;
        fill(h===i&&blink ? 'red' : 'black');
        rect(bx, by+i*s, s, s);
      }
    }
  }
}
