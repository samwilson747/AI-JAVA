import java.util.*;

public class Matrix {

   /**
   @param m    Number of rows.
   @param n    Number of colums.
   @param A    2D Array of doubles holding matrix data.
   */

   private int m, n;
   private int entryCount;
   private double[][] A;


   public Matrix(int rows, int colums) {
      m = rows;
      n = colums;
      entryCount = m * n;
      A = new double[m][n];
      for (int i = 0; i < entryCount; i++) {
        A[i/n][i%n] = 0;
      }
   }

   public Matrix(int rows, int colums, double[] entries) {
      m = rows;
      n = colums;
      entryCount = m * n;
      A = new double[m][n];

      for (int i = 0; i < entryCount; i++) {
         A[i/n][i%n] = entries[i];
      }
   }

   public Matrix(String matrixA) {
     Scanner scan = new Scanner(matrixA);
     if (scan.hasNextInt()) {
         m = scan.nextInt();
     } else {
         System.out.println("Could not read matrix row count.");
         System.exit(1);
     }

     if (scan.hasNextInt()) {
        n = scan.nextInt();
     } else {
        System.out.println("Could not read matrix column count.");
        System.exit(1);
     }

     entryCount = m * n;
     A = new double[m][n];

     for (int i = 0; i < entryCount; i++) {
         if (scan.hasNextDouble()) {
            A[i/n][i%n] = scan.nextDouble();
         } else {
            System.out.println("Could not read matrix entry: " + i);
            System.exit(1);
        }
    }

    scan.close();
   }

   private StringBuilder buildRow(int row) {
      StringBuilder build = new StringBuilder(n*4);

      build.append("| ");

      for (int i = 0; i < n; i++) {
         if (i + 1 == n) {
            build.append(String.format("%.7g", A[row][i]));
         } else {
            build.append(String.format("%.7g", A[row][i]) + ", ");
         }
      }

      build.append(" |\n");
      return build;
   }

   public String debugPrint() {
      StringBuilder build = new StringBuilder(m*n);
      for (int i = 0; i < m; i++){
         build.append(buildRow(i));
      }
      return build.toString();
   }

   public String toString() {
      StringBuilder build = new StringBuilder((entryCount + 1) * 4);
      build.append(m);
      build.append(" ");
      build.append(n);
      build.append(" ");
      for (int i= 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            build.append(String.format("%.2g", A[i][j]));
            if (i + j + 1 < entryCount) {
               build.append(" ");
            }
         }
      }
      return build.toString();
   }
   public int getRowCount() {
      return m;
   }

   public int getColumnCount() {
      return n;
   }

   public double get(int row, int column) {
      return A[row][column];
   }

   public double set(int row, int column, double value) {
      return A[row][column] = value;
   }

   public void output() {
      System.out.print(n + " " + m + " ");
      for(int row = 0; row < m; row++) {
         for (int col = 0; col < n; col++) {
            // if (A[row][col] < 0.000001) {
            //    System.out.print("0.0");
            // }
            // else 
            { System.out.print(A[row][col]); }
            if(row == m -1 && col == n - 1) {
               System.out.println();
            } else {
               System.out.print(" ");
            }
         }
      }
   }

   public static Matrix multiply(Matrix a, Matrix b) {

      int rows = b.getRowCount();
      int columns = a.getColumnCount();

      if (columns != rows) throw new RuntimeException("Illegal matrix dimensions.");
      Matrix c = new Matrix(a.getRowCount(),b.getColumnCount());
      for (int i = 0; i < a.getRowCount(); i++)
          for (int j = 0; j < b.getColumnCount(); j++)
              for (int k = 0; k < a.getColumnCount(); k++)
                  c.set(i, j, (c.get(i,j) + a.get(i, k) * b.get(k, j)));
      return c;
   }

   public static double findColMax(Matrix matrix, int column) {
      double max = 0;
      for(int i = 0; i < matrix.m; i ++) {
         if(matrix.A[i][column] > max)
            max = matrix.A[i][column];
      }
      return max;
   }

   // public static double findRowMax(Matrix matrix, int row) {
   //    double max = 0;
   //    for(int i = 0; i < matrix.m; i ++) {
   //       if(matrix.A[row][i] > max)
   //          max = matrix.A[row][i];
   //    }
   //    return max;
   // }

   public static int findColMaxPos(Matrix matrix, int column) {
      double max = 0;
      int pos = 0;
      for(int i = 0; i < matrix.m; i++) {
         
         if(matrix.A[i][column] - max > 0.000000) {
            max = matrix.A[i][column];
            pos = i;
         }
      }
      return pos;
   }

   public static int[] viterbi(Matrix A, Matrix B, Matrix pi, int[] emSequence) {

      int[] mostProbSeq = new int[emSequence.length];
      Matrix viterbiMatrix = new Matrix(A.m, emSequence.length);
      int[][] viterbiPathMatrix = new int[A.m][emSequence.length];

      // calculate first column in viterbi matrix
      // i == the state being evaluated
      for (int state = 0; state < viterbiMatrix.m; state++) {
         viterbiMatrix.A[state][0] = pi.A[0][state] * B.A[state][emSequence[0]];
         viterbiPathMatrix[state][0] = 0;
      }

      //System.out.println(alphaMatrix.debugPrint());

      for(int e = 1; e < emSequence.length; e++) {
          
         int emission = emSequence[e];
         double maxProb = 0;
         double tempProb = 0;
         int pStateMax = 0;

         // i == the state being evaluated
         for(int state = 0; state < viterbiMatrix.m; state++) {

            maxProb = 0;
            tempProb = 0;
            pStateMax = 0;

            // pState: prior state in veterbi
            for(int pState = 0; pState < viterbiMatrix.m; pState++) {
               
               tempProb = viterbiMatrix.A[pState][e-1] * A.A[pState][state] * B.A[state][emission];
               if (tempProb - maxProb > 0.00000000) {
                  maxProb = tempProb;
                  pStateMax = pState;
               }
            }

            viterbiMatrix.A[state][e] = maxProb;
            viterbiPathMatrix[state][e] = pStateMax;

         }
      }

      System.out.println(viterbiMatrix.debugPrint());

      for(int row = 0; row < A.m; row++) {
         for(int col = 0; col < emSequence.length; col++) {
            System.out.print(viterbiPathMatrix[row][col] + " ");
         }
         System.out.println();
      }
      System.out.println();
      // for(int i = mostProbSeq.length - 1; i >= 0; i--) {
      //    mostProbSeq[i] = findColMaxPos(viterbiMatrix, i);
      // }

      int lastIndex = mostProbSeq.length - 1;
      mostProbSeq[lastIndex] = findColMaxPos(viterbiMatrix, lastIndex);
      //System.out.println(mostProbSeq[viterbiMatrix.n - 1]);

      for(int col = lastIndex - 1; col > 0; col--) {
         mostProbSeq[col] = viterbiPathMatrix[mostProbSeq[col + 1]][col + 1];
      }

      //System.out.println(viterbiPathMatrix[1][3]);

      return mostProbSeq;

   }

   public static void learn(Matrix A, Matrix B, Matrix pi, int[] emSeq) {
      int maxIters = 100;
      int iters = 0;
      double oldLogProb = Double.NEGATIVE_INFINITY;

      Matrix alpha = new Matrix(A.m, emSeq.length);
      double[] alphaScale = new double[emSeq.length];
      Matrix beta  = new Matrix(A.m, emSeq.length);
      double[][] digamma = new double[A.m][emSeq.length];
      double[][][] gamma = new double[A.m][A.m][emSeq.length];

      alphaScale = alphaPass(A, B, pi, emSeq, alpha, alphaScale);
      //System.out.println("Alpha:");
      //System.out.println(alpha.debugPrint());
      betaPass(A, B, emSeq, beta, alphaScale);
      //System.out.println("Beta:");
      //System.out.println(beta.debugPrint());
      gammaCalculations(A, B, emSeq, alpha, beta, gamma, digamma);
      reEstimate(A, B, pi, emSeq, gamma, digamma);
      double logProb = logProb(alphaScale);

      // System.out.println("Made it through the computations");
      // System.out.println("A:");
      // System.out.println(A.debugPrint());
      // System.out.println("B:");
      // System.out.println(B.debugPrint());

      iters = iters + 1;
      
      while (iters < maxIters && logProb > oldLogProb) {
         oldLogProb = logProb;
         // not done learning... run again.
         alphaScale = alphaPass(A, B, pi, emSeq, alpha, alphaScale);
         betaPass(A, B, emSeq, beta, alphaScale);
         gammaCalculations(A, B, emSeq, alpha, beta, gamma, digamma);
         reEstimate(A, B, pi, emSeq, gamma, digamma);
         logProb = logProb(alphaScale);

         iters++;
         // System.out.println("Learning Pass: " + (iters + 1) + " Finished");
         // System.out.println("OldLogProb: " + oldLogProb);
         // System.out.println("NewLogProb: " + logProb);
         
      }

      // output learned A and B 
      // System.out.println("Final Results");
      // System.out.println("A:");
      // System.out.println(A.debugPrint());
      // System.out.println("B:");
      // System.out.println(B.debugPrint());

      A.output();
      B.output();
   }

   public static double[] alphaPass(Matrix A, Matrix B, Matrix pi, int[] emSeq, Matrix alpha, double[] scale) {

      //compute first alpha column
      for(int state = 0; state < A.m; state++) {
         alpha.A[state][0] = pi.A[0][state]*B.A[state][emSeq[0]];
         scale[0] = scale[0] + alpha.A[state][0];
      }

      //scale  first alpha column
      scale[0] = 1.0 / scale[0];
      for(int state = 0; state < A.m; state++) {
         alpha.A[state][0] = scale[0] * alpha.A[state][0];
      }

      //compute alpha
      for(int time = 1; time < emSeq.length; time++) {

         scale[time] = 0;

         for(int state = 0; state < A.m; state++) {

            alpha.A[state][time] = 0;

            for(int pState = 0; pState < A.m; pState++) {
               alpha.A[state][time] = alpha.A[state][time] + alpha.A[pState][time - 1] * A.A[pState][state];
            }

            alpha.A[state][time] = alpha.A[state][time] * B.A[state][emSeq[time]];
            scale[time] = scale[time] + alpha.A[state][time];

         }

         //scale alpha at time t
         scale[time] = 1.0 / scale[time];
         for(int state = 0; state < A.m; state++) {
            alpha.A[state][time] = scale[time] * alpha.A[state][time];
         }
      }

      return scale;
   }

   public static void betaPass(Matrix A, Matrix B, int[] emSeq, Matrix beta, double[] alphaScale) {

      // Let beta[state][time - 1] = 1, scaled by scale[time - 1]
      for(int state = 0; state < beta.m; state++) {
         beta.A[state][beta.n - 1] = alphaScale[beta.n - 1];
      }

      // System.out.println("alphaScale:");
      // for (double val : alphaScale) {
      //    System.out.print(val+" ");
      // }
      // System.out.println("\nbeta:");
      // System.out.println(beta.debugPrint());

      // beta-pass
      for(int time = alphaScale.length - 2; time >= 0; time--) {

         for(int state = 0; state < A.m; state++) {

            beta.A[state][time] = 0;
            for(int nState = 0; nState < A.m; nState++) {
               double temp = A.A[state][nState] * B.A[nState][emSeq[time + 1]] * beta.A[nState][time + 1];
               beta.A[state][time] = beta.A[state][time] + temp;
            }

            // scale beta[state][time] with same scale factor as alpha[state][time]
            beta.A[state][time] = alphaScale[time]*beta.A[state][time];
         }
      }
   }

   public static void gammaCalculations(Matrix A, Matrix B, int[] emSeq, Matrix alpha, Matrix beta, double[][][] gamma, double[][] digamma) {

      for(int t = 0; t < emSeq.length - 1; t++) {
         double temp = 0;
         double denom = 0;
         for(int state = 0; state < A.m; state++) {
            digamma[state][t] = 0;
            for(int nState = 0; nState < A.n; nState++) {
               temp = alpha.A[state][t] * A.A[state][nState] * B.A[nState][emSeq[t + 1]] * beta.A[nState][t + 1];
               denom = denom + temp;
            }
         }

         for(int state = 0; state < A.m; state++) {

            //digamma[state][t] = 0;
            for(int nState = 0; nState < A.n; nState++) {
               temp = alpha.A[state][t] * A.A[state][nState] * B.A[nState][emSeq[t + 1]] * beta.A[nState][t + 1];
               gamma[state][nState][t] = temp/denom;
               digamma[state][t] = digamma[state][t] + gamma[state][nState][t];
            }
         }
      }

      //special case for diGamma[state][time - 1]
      double denom = 0;
      int lastIndex = emSeq.length - 1;

      for(int state = 0; state < A.n; state++) {
         denom = denom + alpha.A[state][lastIndex];
      }

      for(int state = 0; state < A.n; state++) {
         digamma[state][lastIndex] = alpha.A[state][lastIndex]/denom;
      }
   }

   public static void reEstimate(Matrix A, Matrix B, Matrix pi, int[] emSeq, double[][][] gamma, double[][] digamma) {
      //re-estimate pi
      for(int state = 0; state < pi.n; state++) {
         pi.A[0][state] = digamma[state][0];
      }

      //re-estimate A
      for (int state = 0; state < A.n; state++) {
         for (int nState = 0; nState < A.n; nState++) {

            double numer = 0;
            double denom = 0;
            //maybe needs to be emSeq.length - 1
            for (int t = 0; t < emSeq.length - 1; t++) {
               numer = numer + gamma[state][nState][t];
               denom = denom + digamma[state][t];
            }
            
            A.A[state][nState] = numer/denom;
         }
      }

      //re-estimate B
      for (int state = 0; state < B.m; state++) {
         for(int emission = 0; emission < B.n; emission++) {
            double numer = 0;
            double denom = 0;

            for(int t = 0; t < emSeq.length; t++) {
               if (emSeq[t] == emission) {
                  numer = numer + digamma[state][t];
               }
               denom = denom + digamma[state][t];
            }

            B.A[state][emission] = numer/denom;
         }
      }
   }

   public static double logProb(double[] scale) {

      double logProb = 0;
      for (int time = 0; time < scale.length; time++) {
         logProb = logProb + Math.log(scale[time]);
      }

      return -logProb;
   }

}