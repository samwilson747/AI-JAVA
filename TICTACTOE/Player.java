import java.util.*;

public class Player {
   /**
   * Performs a move
   *
   * @param gameState
   *            the current state of the board
   * @param deadline
   *            time before which we must have returned
   * @return the next state the board is in after our move
   */

   public int evaluation(GameState state, int player, int depth) {
      if(state.isXWin())
         return 10;//*(depth+1);
      else if(state.isOWin())
         return -10;//*(depth+1);
      else{
         return 0;
         /*
         int count = 0;
         Vector<GameState> children = new Vector<GameState>();
         state.findPossibleMoves(children);
         for(int i = 0; i < children.size(); i++)
         {
            if( children.elementAt(i).isXWin())
               count++;
         }
         return count*depth; */
      }

   }

   /*
   public int minimax(GameState state, int depth, int player) {
      Vector<GameState> children = new Vector<GameState>();
      state.findPossibleMoves(children);

      if(depth == 0 || state.isEOG()) {
         if(state.isXWin())
            return 10;
         else if(state.isOWin())
            return -10;
         else
            return -5;
      }
      else {
         if(player == Constants.CELL_X) {
            int bp = -10;
            for(int i = 0;i<children.size();i++) {
               int val = minimax( children.elementAt(i) , depth-1, children.elementAt(i).getNextPlayer() );
               bp = Math.max(bp, val);
            }
            return bp;
         }
         else {
            int bp = 10;
            for(int i = 0;i<children.size();i++) {
               int val = minimax(children.elementAt(i), depth-1, children.elementAt(i).getNextPlayer());
               bp = Math.max(bp, val);
            }
            return bp;
         }
      }
   }
   */
   public int alphabeta(GameState state, int depth, int alpha, int beta, int player) {
      Vector<GameState> children = new Vector<GameState>();
      state.findPossibleMoves(children);
      int v;

      if( depth == 0 || state.isEOG()) {
         return evaluation(state, player, depth);
      }
      else if(player == Constants.CELL_O) { //PLAYER MAX
         v = Integer.MIN_VALUE;
         for(int i = 0; i < children.size(); i++) {
            GameState currChild = children.elementAt(i);
            v = Math.max(v, alphabeta(currChild, depth-1, alpha, beta, state.getNextPlayer()));
            alpha = Math.max(alpha, v);
            if(beta <= alpha)
               break; //tree pruned
         }
      }
      else { //PLAYER MIN
         v = Integer.MAX_VALUE;
         for(int i = 0; i < children.size(); i++) {
            GameState currChild = children.elementAt(i);
            v = Math.min(v, alphabeta(currChild, depth-1, alpha, beta, state.getNextPlayer()));
            beta = Math.min(beta, v);
            if(beta<= alpha)
               break; //tree pruned
         }
      }
      return v;
   }
   public GameState play(final GameState gameState, final Deadline deadline) {
      Vector<GameState> nextStates = new Vector<GameState>();
      gameState.findPossibleMoves(nextStates);

      if (nextStates.size() == 0) {
         // Must play "pass" move if there are no other moves possible.
         return new GameState(gameState, new Move());
      }
      int idx, max;
      idx = 0;
      max = -10;
      int alpha = Integer.MIN_VALUE;
      int beta = Integer.MAX_VALUE;
      int depth = 5;

      for(int i = 0; i < nextStates.size(); i++) {
         //int val = minimax(nextStates.elementAt(i), depth, Constants.CELL_X);
         int val = alphabeta(nextStates.elementAt(i), depth, alpha, beta, gameState.getNextPlayer());
         if(val > max)
            idx = i;
      }
      return nextStates.elementAt(idx);
      /**
      * Here you should write your algorithms to get the best next move, i.e.
      * the best next state. This skeleton returns a random move instead.
      */
   //   Random random = new Random();
   //   return nextStates.elementAt(random.nextInt(nextStates.size()));
   }
}
