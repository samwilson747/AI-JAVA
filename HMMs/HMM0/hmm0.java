import java.util.Scanner;

public class hmm0 {
   public static void main(String[] args) {

      Matrix tMatrix = null;
      Matrix eMatrix = null;
      Matrix ispMatrix = null;

      Scanner scan = new Scanner(System.in);

      // read transition matrix
      if (scan.hasNextLine()) {
         String tMatrixString = scan.nextLine();
         tMatrix = new Matrix(tMatrixString);
      } else {
         //System.out.println("No Transition Matrix Found!");
         System.exit(1);
      }

      // read emission matrix
      if (scan.hasNextLine()) {
         String eMatrixString = scan.nextLine();
         eMatrix = new Matrix(eMatrixString);
      } else {
         //System.out.println("No Emission Matrix Found!");
         System.exit(1);
      }

      // read initial state probability matrix

      if (scan.hasNextLine()) {
         String ispMatrixString = scan.nextLine();
         ispMatrix = new Matrix(ispMatrixString);
      } else {
         //System.out.println("No Initial State Probability Matrix Found!");
         System.exit(1);
      }

      scan.close();

      Matrix result = Matrix.multiply(ispMatrix, tMatrix);
      Matrix output = Matrix.multiply(result, eMatrix);

      System.out.println(output);
   }
}
