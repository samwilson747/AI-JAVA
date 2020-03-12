double alphaPass(const Matrix &A, const Matrix &B, const vector<double> &pi,
 const vector<int> &seq) {
   int obs, state, nxtState;
   const int kStates = pi.size(), kObservations = seq.size();
   double probSum, ans;
   Matrix alpha = Matrix(kObservations, kStates); //Observations * states

   for (state = 0; state < kStates; state++)
      alpha[0][state] = B[state][seq[0]] * pi[state];

   for (obs = 1; obs < kObservations; obs++){
      for (state = 0; state < kStates; state++){
         probSum = 0.0;
         for (nxtState = 0; nxtState < kStates; nxtState++){
            probSum += A[nxtState][state] * alpha[obs-1][nxtState];
         }
         alpha[obs][state] = probSum * B[state][seq[obs]];
      }
   }

   ans = 0.0;
   for (state = 0; state < kStates; state++)
      ans += alpha[kObservations-1][state];
   return ans;
}
