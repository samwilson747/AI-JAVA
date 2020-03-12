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
      double max = Double.NEGATIVE_INFINITY;
      int pos = -1;
      for(int i = 0; i < matrix.m; i++) {
         
         if(matrix.A[i][column] > max) {
            max = matrix.A[i][column];
            pos = i;
         }
      }
      return pos;
   }

   public static String viterbi(Matrix A, Matrix B, Matrix pi, int[] emissionSequence) {

      Matrix viterbiMatrix = new Matrix(A.m, emissionSequence.length);
      int[][] viterbiPathMatrix = new int[A.m][emissionSequence.length];
      int lastEmissionIndex = emissionSequence.length - 1;
      double lastStateMaxProb = Double.NEGATIVE_INFINITY;
      int lastStateMaxIndex = -1;

      // calculate first column in viterbi matrix
      // i == the state being evaluated
      for (int state = 0; state < viterbiMatrix.m; state++) {
         viterbiMatrix.A[state][0] = Math.log(pi.A[0][state] * B.A[state][emissionSequence[0]]);
         viterbiPathMatrix[state][0] = 0;
      }

      //System.out.println(alphaMatrix.debugPrint());

      for(int e = 1; e < emissionSequence.length; e++) {
          
         int emission = emissionSequence[e];
         double maxProb = 0;
         double tempProb = 0;
         int pStateMax = 0;

         // i == the state being evaluated
         for(int state = 0; state < viterbiMatrix.m; state++) {

            viterbiMatrix.A[state][e] = Double.NEGATIVE_INFINITY;
            tempProb = 0;

            // pState: prior state in veterbi
            for(int pState = 0; pState < viterbiMatrix.m; pState++) {
               
               tempProb = viterbiMatrix.A[pState][e-1] + Math.log(A.A[pState][state]) + Math.log(B.A[state][emission]);
               if (tempProb > viterbiMatrix.A[state][e]) {
                  viterbiMatrix.A[state][e] = tempProb;
                  viterbiPathMatrix[state][e] = pState;
               }
            }

            if (e == lastEmissionIndex && lastStateMaxProb < viterbiMatrix.A[state][e]) {
               lastStateMaxProb = viterbiMatrix.A[state][e];
               lastStateMaxIndex = state;
            }
         }
      }

      // System.out.println(viterbiMatrix.debugPrint());

      // for(int row = 0; row < A.m; row++) {
      //    for(int col = 0; col < emissionSequence.length; col++) {
      //       System.out.print(viterbiPathMatrix[row][col] + " ");
      //    }
      //    System.out.println();
      // }


      String finalPath = lastStateMaxIndex + "";

      for(int t = lastEmissionIndex; t > 0; t--) {

         lastStateMaxIndex = viterbiPathMatrix[lastStateMaxIndex][t];
         finalPath = lastStateMaxIndex + " " + finalPath;
      }

      //System.out.println(viterbiPathMatrix[1][3]);

      return finalPath;

   }
}