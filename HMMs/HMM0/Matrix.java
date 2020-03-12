import java.util.Scanner;

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
            build.append(A[row][i]);
         } else {
            build.append(A[row][i] + ", ");
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
}
