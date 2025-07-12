import controlP5.*;
import java.io.File;

ControlP5 cp5;
int numElements = 5;
float barLength = 1.0f;
float E = 210e9f;
float A = 0.01f;

float[] nodeX;
float[] displacement;
float[][] K;
float[] F;
float[] stress;
float[] strain;

int highlightElement = -1;
int highlightNode = -1;

void settings() {
  size(1000, 700);
}

void setup() {
  cp5 = new ControlP5(this);

  cp5.addSlider("numElements")
     .setPosition(20, 20)
     .setRange(2, 20)
     .setValue(numElements)
     .setNumberOfTickMarks(19)
     .setSliderMode(Slider.FLEXIBLE);

  cp5.addTextfield("highlightElement")
     .setPosition(20, 60)
     .setSize(100, 30)
     .setLabel("Enter Node #")
     .setAutoClear(false)
     .onChange(e -> {
       try {
         highlightNode = int(e.getController().getStringValue());
         println("Showing connectivity for Node: " + highlightNode);
       } catch (Exception ex) {
         highlightNode = -1;
       }
     });

  computeFEM();
}

void draw() {
  background(255);
  drawStiffnessMatrix();
  drawDisplacementPlot();
  fill(0);
  text("Press 'E' to export results", 20, height - 20);
  text("Node input highlights elements connected to that node", 20, height - 40);
}

void keyPressed() {
  if (key == 'e' || key == 'E') {
    exportResults();
  }
}

void computeFEM() {
  int NN = numElements + 1;
  nodeX = new float[NN];
  displacement = new float[NN];
  F = new float[NN];
  K = new float[NN][NN];
  stress = new float[numElements];
  strain = new float[numElements];

  float dx = barLength / numElements;
  for (int i = 0; i < NN; i++) {
    nodeX[i] = i * dx;
    F[i] = 0;
  }
  F[NN - 1] = 1000.0f;

  for (int e = 0; e < numElements; e++) {
    int i = e;
    int j = e + 1;
    float Le = nodeX[j] - nodeX[i];
    float k = E * A / Le;

    K[i][i] += k;
    K[j][j] += k;
    K[i][j] -= k;
    K[j][i] -= k;
  }

  float big = 1e20f;
  K[0][0] += big;
  F[0] += big * 0;

  displacement = solveLinearSystem(K, F);

  for (int e = 0; e < numElements; e++) {
    int i = e;
    int j = e + 1;
    float Le = nodeX[j] - nodeX[i];
    strain[e] = (displacement[j] - displacement[i]) / Le;
    stress[e] = E * strain[e];
  }

  printStiffnessMatrix();
}

void printStiffnessMatrix() {
  println("==== Global Stiffness Matrix ====");
  for (int i = 0; i < K.length; i++) {
    String row = "";
    for (int j = 0; j < K[i].length; j++) {
      row += nf(K[i][j], 0, 2) + "\t";
    }
    println(row);
  }
}

float[] solveLinearSystem(float[][] A, float[] b) {
  int n = b.length;
  float[][] M = new float[n][n];
  float[] x = new float[n];
  float[] f = new float[n];

  for (int i = 0; i < n; i++) {
    f[i] = b[i];
    for (int j = 0; j < n; j++) M[i][j] = A[i][j];
  }

  for (int k = 0; k < n - 1; k++) {
    for (int i = k + 1; i < n; i++) {
      float factor = M[i][k] / M[k][k];
      for (int j = k; j < n; j++) M[i][j] -= factor * M[k][j];
      f[i] -= factor * f[k];
    }
  }

  for (int i = n - 1; i >= 0; i--) {
    float sum = 0;
    for (int j = i + 1; j < n; j++) sum += M[i][j] * x[j];
    x[i] = (f[i] - sum) / M[i][i];
  }
  return x;
}

void drawStiffnessMatrix() {
  int s = 10;
  int NN = numElements + 1;
  pushMatrix();
  translate(250, 20);
  for (int i = 0; i < NN; i++) {
    for (int j = 0; j < NN; j++) {
      float val = K[i][j];
      if (abs(val) > 1e-6f) {
        if (highlightNode >= 0 && (i == highlightNode || j == highlightNode)) {
          fill(255, 0, 0);
        } else {
          fill(0);
        }
        rect(j * s, i * s, s, s);
      }
    }
  }
  popMatrix();
  fill(0);
  text("Stiffness Matrix (Red = connected to input node)", 250, NN * 10 + 40);
}

void drawDisplacementPlot() {
  float x0 = 250;
  float y0 = height - 100;
  float scaleX = 500;
  float scaleY = 1e5f;

  stroke(0);
  fill(0);
  text("Displacement Plot", 20, height - 120);
  noFill();
  beginShape();
  for (int i = 0; i < nodeX.length; i++) {
    float x = x0 + nodeX[i] * scaleX;
    float y = y0 - displacement[i] * scaleY;
    vertex(x, y);
  }
  endShape();
}

void exportResults() {
  String folderPath = sketchPath("fem_output");
  File folder = new File(folderPath);
  if (!folder.exists()) folder.mkdirs();

  PrintWriter writer = createWriter(folderPath + "/fem1d_results.txt");

  writer.println("==== Node Displacements (in meters) ====");
  for (int i = 0; i < displacement.length; i++) {
    writer.println("Node " + i + ": " + displacement[i]);
  }

  writer.println("\n==== Element Strain ====");
  for (int i = 0; i < strain.length; i++) {
    writer.println("Element " + i + ": " + strain[i]);
  }

  writer.println("\n==== Element Stress (Pascals) ====");
  for (int i = 0; i < stress.length; i++) {
    writer.println("Element " + i + ": " + stress[i]);
  }

  writer.println("\n==== Global Stiffness Matrix ====");
  for (int i = 0; i < K.length; i++) {
    for (int j = 0; j < K[i].length; j++) {
      writer.print(nf(K[i][j], 0, 2) + "\t");
    }
    writer.println();
  }

  writer.flush();
  writer.close();
  println("✅ Exported all results to: fem_output/fem1d_results.txt");
}

void controlEvent(ControlEvent theEvent) {
  if (theEvent.isFrom("numElements")) {
    numElements = int(theEvent.getController().getValue());
    computeFEM();
  }
}
EEEE
