import java.util.Scanner;

public class hmm2 {
   public static void main(String[] args) {

      Matrix A = null;
      Matrix B = null;
      Matrix pi = null;
      Matrix alpha = null;
      int[] sequence = null;

      Scanner scan = new Scanner(System.in);

      // read transition matrix
      if (scan.hasNextLine()) {
         String AString = scan.nextLine();
         A = new Matrix(AString);
      } else {
         //System.out.println("No Transition Matrix Found!");
         System.exit(1);
      }

      // read emission matrix
      if (scan.hasNextLine()) {
         String BString = scan.nextLine();
         B = new Matrix(BString);
      } else {
         //System.out.println("No Emission Matrix Found!");
         System.exit(1);
      }

      // read initial state probability matrix

      if (scan.hasNextLine()) {
         String piString = scan.nextLine();
         pi = new Matrix(piString);
      } else {
         //System.out.println("No Initial State Probability Matrix Found!");
         System.exit(1);
      }

      // read in the emission sequence
      if (scan.hasNextLine()) {
         String emissionSeqString = scan.nextLine();
         Scanner emScan = new Scanner(emissionSeqString);
         if (emScan.hasNextInt()) {
            sequence = new int[emScan.nextInt()];
            for (int i = 0; i < sequence.length; i++) {
               if (emScan.hasNextInt()) {
                  sequence[i] = emScan.nextInt();
               } else {
                  //System.out.println("Could not read emission sequence at value: " + i);
                  System.exit(1);
               }
            }

         } else {
            //System.out.println("No Emission Sequence Found!!!!");
            System.exit(1);
         }

      } else {
         //System.out.println("No Emission Sequence Found!");
         System.exit(1);
      }

      scan.close();

      System.out.println(Matrix.viterbi(A, B, pi, sequence));
      
   }
   
}
