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

   // public static int findColMaxPos(Matrix matrix, int column) {
   //    double max = 0;
   //    int pos = 0;
   //    for(int i = 0; i < matrix.m; i ++) {
   //       if(matrix.A[i][column] > max)
   //          max = matrix.A[i][column];
   //          pos = i;
   //    }
   //    return pos;
   // }

   public double findTotalEmissionSeqProbability() {
      double max = 0;
      for (int i = 0; i < m; i++) {
         max += A[i][n - 1];
      }

      return max;
   }

   public static Matrix computeAlphaMatrix(Matrix A, Matrix B, Matrix pi, int[] emisionSequence) {

      double totalProb = 0;
      Matrix alphaMatrix = new Matrix(A.m, emisionSequence.length);

      // calculate first column in alpha matrix
      // i == the state being evaluated
      for (int i = 0; i < alphaMatrix.m; i++) {
          alphaMatrix.A[i][0] = pi.A[0][i] * B.A[i][emisionSequence[0]];
          //alphaMatrix.A[i][0] *= B.A[i][emisionSequence[0]];
      }

      //System.out.println(alphaMatrix.debugPrint());

      for(int e = 1; e < emisionSequence.length; e++) {
          
         int emission = emisionSequence[e];

         // i == the state being evaluated
         for(int i = 0; i < alphaMatrix.m; i++) {

            totalProb = 0;

            for(int j = 0; j < alphaMatrix.m; j++) {
               totalProb += alphaMatrix.A[j][e-1] * A.A[j][i];
            }

            alphaMatrix.A[i][e] = totalProb * B.A[i][emission];
        }
      }

      //System.out.println(alphaMatrix.debugPrint());
      return alphaMatrix;

   }
}