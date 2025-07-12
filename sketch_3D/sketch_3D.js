// === 3D Hexahedral FEM VISUALIZER ===
// 8‑node trilinear hexahedral elements, frontal solver stub

let nx=1, ny=1, nz=1;
let fem;
let sliderX, sliderY, sliderZ, dofInput, showButton;
let highlightDOF=-1;

function setup(){
  createCanvas(1400,800,WEBGL);
  textFont('monospace');
  sliderX = createSlider(1,3,nx,1).position(10,10);
  sliderY = createSlider(1,3,ny,1).position(10,40);
  sliderZ = createSlider(1,3,nz,1).position(10,70);
  dofInput = createInput('').position(10,100).size(60);
  showButton = createButton('🔎 Highlight DOF').position(80,100)
    .mousePressed(()=>{ highlightDOF=int(dofInput.value()); });
  buildFEM();
}

function draw(){
  background(255);
  orbitControl();
  fem.drawMesh();
  push();
    resetMatrix(); ortho(); translate(-width/2,-height/2);
    fem.drawMatrix(highlightDOF);
    fem.drawForceVector(highlightDOF);
  pop();
  fill(0);
  text(`nx ${nx}, ny ${ny}, nz ${nz}`,10,140);
}

function buildFEM(){
  nx=sliderX.value(); ny=sliderY.value(); nz=sliderZ.value();
  fem = new HexaFEM(nx,ny,nz);
  fem.assemble();
}

class HexaFEM{
  constructor(nx,ny,nz){
    this.nx=nx; this.ny=ny; this.nz=nz;
    this.nodes=[]; this.elements=[];
    this.buildGrid();
    this.DOF=this.nodes.length*3;
    this.K=Array(this.DOF).fill().map(()=>Array(this.DOF).fill(0));
    this.F=Array(this.DOF).fill(0);
    this.E=2e5; this.nu=0.3;
  }
  buildGrid(){
    for(let k=0;k<=this.nz;k++){
      for(let j=0;j<=this.ny;j++){
        for(let i=0;i<=this.nx;i++){
          this.nodes.push([i,j,k]);
        }} }
    let sx=this.nx+1, sy=this.ny+1;
    for(let k=0;k<this.nz;k++)for(let j=0;j<this.ny;j++)for(let i=0;i<this.nx;i++){
      let n0 = k*sx*sy + j*sx + i;
      let n1=n0+1, n2=n0+sx+1, n3=n0+sx;
      let offset = sx*sy;
      let n4=n0+offset, n5=n1+offset, n6=n2+offset, n7=n3+offset;
      this.elements.push([n0,n1,n2,n3,n4,n5,n6,n7]);
    }
  }
  assemble(){
    // stub K,F fill
    this.elements.forEach(el=>{
      el.forEach(n=>{
        let idx = 3*n;
        this.F[idx+2] = -80000; // apply Z-load
      });
    });
    // penalty BC on node 0
    [0].forEach(n=>{ for(let d=0;d<3;d++) this.K[3*n+d][3*n+d]+=1e20; });
  }
  drawMesh(){
    stroke(0); noFill();
    this.elements.forEach(el=>{
      beginShape();
      el.forEach(n=>{
        let [i,j,k]=this.nodes[n];
        vertex(i*50,j*50,k*50);
      }); endShape(CLOSE);
    });
    stroke(255,0,0); this.nodes.forEach((p,i)=>{
      push(); translate(p[0]*50,p[1]*50,p[2]*50);
      sphere(3); pop(); });
  }
  drawMatrix(h){
    let s=6, bx=600, by=50;
    for(let i=0;i<this.DOF;i++)for(let j=0;j<this.DOF;j++){
      if(abs(this.K[i][j])>1e-3){
        let blink=frameCount%60<30;
        fill((h===i||h===j)&&blink?'red':'black');
        rect(bx+j*s,by+i*s,s,s);
      }
    }
  }
  drawForceVector(h){
    let s=6, bx=600+this.DOF*6+20, by=50;
    for(let i=0;i<this.DOF;i++){
      if(abs(this.F[i])>1e-3){
        let blink=frameCount%60<30;
        fill(h===i&&blink?'red':'black');
        rect(bx,by+i*s,s,s);
      }
    }
  }
}
