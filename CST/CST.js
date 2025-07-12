// === 2D CST FEM VISUALIZER WITH DYNAMIC MESH GENERATION ===
// Implements Constant Strain Triangle elements with automatic grid mesh

let nx = 2, ny = 2;  // number of quads in x and y
let fem;
let sliderX, sliderY, dofInput, showButton;
let highlightDOF = -1;

function setup() {
  createCanvas(1280, 720);
  textFont('monospace');

  sliderX = createSlider(1, 10, nx, 1);
  sliderX.position(20, 20);
  sliderX.input(resetFEM);
  sliderY = createSlider(1, 10, ny, 1);
  sliderY.position(20, 50);
  sliderY.input(resetFEM);

  dofInput = createInput('');
  dofInput.position(20, 90).size(60);
  showButton = createButton('🔎 Highlight DOF');
  showButton.position(90, 90).mousePressed(() => highlightDOF = int(dofInput.value()));

  resetFEM();
}

function draw() {
  background(255);
  fem.drawMesh();
  fem.drawMatrix(highlightDOF);
  fem.drawForceVector(highlightDOF);
  fem.drawElementStiffness();

  fill(0);
  text(`nx: ${sliderX.value()}  ny: ${sliderY.value()}`, 20, 140);
}

function resetFEM() {
  nx = sliderX.value(); ny = sliderY.value();
  fem = new CSTFEM(nx, ny);
  fem.assemble();
}

class CSTFEM {
  constructor(nx, ny) {
    this.nx = nx; this.ny = ny;
    this.nodes = [];
    this.elements = [];
    this.buildMesh();
    this.DOF = this.nodes.length * 2;
    this.K = Array(this.DOF).fill().map(() => Array(this.DOF).fill(0));
    this.F = Array(this.DOF).fill(0);
    this.E = 3e7; this.nu = 0.25; this.t = 1.0;
    this.SEcollection = [];
  }

  buildMesh() {
    let w = 400, h = 400;
    for (let j = 0; j <= this.ny; j++) {
      for (let i = 0; i <= this.nx; i++) {
        let x = 200 + i * (w/this.nx);
        let y = 200 + j * (h/this.ny);
        this.nodes.push([x,y]);
      }
    }
    for (let j = 0; j < this.ny; j++) {
      for (let i = 0; i < this.nx; i++) {
        let n0 = j*(this.nx+1)+i;
        let n1 = n0+1;
        let n2 = n0+this.nx+1;
        let n3 = n2+1;
        // two triangles
        this.elements.push([n0,n1,n3]);
        this.elements.push([n0,n3,n2]);
      }
    }
  }

  assemble() {
    const D = this.getDMatrix();
    for (let e=0; e<this.elements.length; e++) {
      let [i1,i2,i3]=this.elements[e];
      let p1=this.nodes[i1], p2=this.nodes[i2], p3=this.nodes[i3];
      let A = 0.5*abs((p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]));
      let B = this.getBMatrix(p1,p2,p3,A);
      let SE = Array(6).fill().map(()=>Array(6).fill(0));
      for (let i=0;i<6;i++) for(let j=0;j<6;j++) for(let k=0;k<3;k++) for(let l=0;l<3;l++)
        SE[i][j]+=B[k][i]*D[k][l]*B[l][j]*A*this.t;
      let dof=[2*i1,2*i1+1,2*i2,2*i2+1,2*i3,2*i3+1];
      this.SEcollection.push({e,dof,SE});
      for (let i=0;i<6;i++) for(let j=0;j<6;j++) this.K[dof[i]][dof[j]]+=SE[i][j];
    }
    // apply BC at bottom-left node
    this.K[0][0]+=1e20; this.K[1][1]+=1e20;
    // apply loads at top-right node
    let n= this.nodes.length-1;
    this.F[2*n]+= -200; this.F[2*n+1]+= -400;
  }

  getDMatrix() {
    let E=this.E, nu=this.nu;
    let c=E/(1-nu*nu);
    return [[c, c*nu,0],[c*nu,c,0],[0,0,c*(1-nu)/2]];
  }

  getBMatrix(p1,p2,p3,A) {
    let [x1,y1]=p1,[x2,y2]=p2,[x3,y3]=p3;
    let b1=y2-y3,b2=y3-y1,b3=y1-y2;
    let c1=x3-x2,c2=x1-x3,c3=x2-x1;
    let B=Array(3).fill().map(()=>Array(6).fill(0));
    B[0][0]=b1;B[0][2]=b2;B[0][4]=b3;
    B[1][1]=c1;B[1][3]=c2;B[1][5]=c3;
    B[2][0]=c1;B[2][1]=b1;B[2][2]=c2;B[2][3]=b2;B[2][4]=c3;B[2][5]=b3;
    for(let i=0;i<3;i++)for(let j=0;j<6;j++)B[i][j]/=(2*A);
    return B;
  }

  drawMesh() {
    stroke(200); fill(0);
    for(let i=0;i<this.nodes.length;i++){
      let [x,y]=this.nodes[i]; text(i,x+5,y-5);
      ellipse(x,y,4);
    }
    for(let e of this.elements){let p1=this.nodes[e[0]],p2=this.nodes[e[1]],p3=this.nodes[e[2]];
      noFill(); stroke(150); triangle(p1[0],p1[1],p2[0],p2[1],p3[0],p3[1]);
    }
  }

  drawMatrix(h) {
    let s=10, bx=600, by=50;
    for(let i=0;i<this.DOF;i++)for(let j=0;j<this.DOF;j++){
      if(abs(this.K[i][j])>1e-3){fill((h===i||h===j)&&frameCount%60<30?255:0);
        rect(bx+j*s,by+i*s,s,s);} }
  }

  drawForceVector(h) {
    let s=10, bx=600+this.DOF*10+20, by=50;
    for(let i=0;i<this.DOF;i++)if(abs(this.F[i])>1e-3){fill((h===i)&&frameCount%60<30?255:0);
      rect(bx,by+i*10,10,10);} }

  drawElementStiffness() {
    let bx=20, by=360;
    for(let k=0;k<this.SEcollection.length;k++){
      let {e,SE}=this.SEcollection[k]; fill(0);
      text(`E${e}`,bx+k*140,by);
      for(let i=0;i<6;i++)for(let j=0;j<6;j++){
        text(SE[i][j].toExponential(2),bx+k*140+j*20,by+15+i*12);
      }
    }
  }
}
