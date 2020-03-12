/*hmm0.cpp
 *
 *
 *
 * Vectors were used instead of arrays because of:
 * https://stackoverflow.com/questions/10865861/when-to-use-vectors-and-when-to-use-arrays-in-c
 */

#include <iostream>
#include <vector>
#include <sstream>
#include <iterator>
#include <limits>
#include <math.h>

//-----------Global Variables--------------
const int LOGZERO = 1;
int states = -1;//Set in the function readInput
int measurments = -1;//Set in the function readInput
int tMax = -1;//Set in the function readInput

//----------Function Declarations----------
void readInput(std::vector<std::vector<double>>* A,std::vector<std::vector<double>>* B, std::vector<double>* pi, std::vector<double>* seq);
double alphaLDPass(std::vector<std::vector<double>>* A,std::vector<std::vector<double>>* B, std::vector<double>* pi, std::vector<double>* seq, std::vector<std::vector<double>>* alphaLD);
void betaLDPass(std::vector<std::vector<double>>* A, std::vector<std::vector<double>>* B, std::vector<double>* seq, std::vector<std::vector<double>>* beta);
double getGammaLDDenominator(std::vector<std::vector<double>>* alpha);
void gammaLDFun(std::vector<std::vector<double>>* alpha, std::vector<std::vector<double>>* beta, std::vector<std::vector<double>>* gamma);
void diGammaLDFun(std::vector<std::vector<double>>* A, std::vector<std::vector<double>>* B, std::vector<double>* seq, std::vector<std::vector<double>>* alpha, std::vector<std::vector<double>>* beta, std::vector<std::vector<std::vector<double>>>* diGamma);
void reEstimateParamsLD(std::vector<std::vector<double>>* A, std::vector<std::vector<double>>* B, std::vector<double>* pi, std::vector<double>* seq, std::vector<std::vector<double>>* gamma, std::vector<std::vector<std::vector<double>>>* diGamma);
void outputMatrix(std::vector<std::vector<double>> matrix);
double eexp(double x);
double eln(double x);
double elnsum(double elnx, double elny);
double elnproduct(double elnx,double elny);
void displayMatrixLD(std::vector<std::vector<double>> matrix);
void displayVector(std::vector<double> vector);

//----------Function Definitions-----------

//-Core Functions
//-Log Helper Functions
//-Debugging Functions
//-Main Function

//--------Core Functions-------
void readInput(std::vector<std::vector<double>>* A,std::vector<std::vector<double>>* B, std::vector<double>* pi, std::vector<double>* seq){
     //Reads the input and store it in the transition matrix, emission matrix, initialization vector and a sequence vector
     int matrix = 0;//0=A, 1=B, 2=pi

     for (std::string line; std::getline(std::cin, line);) {
         std::vector<double> numbers;

         //Convert line to vector of doubles
         std::istringstream iss(line);
         std::copy(std::istream_iterator<double>(iss),
                   std::istream_iterator<double>(),
                   std::back_inserter(numbers));

         switch(matrix){
             case 0://Resize transition matrix and fill in data
             {
                 states = numbers[0];
                 numbers.erase(numbers.begin());
                 numbers.erase(numbers.begin());

                 (*A).resize(states);
                 for(int i = 0; i < states; i++){
                     (*A)[i].resize(states);
                 }

                 for(int row = 0; row < states; row++){
                     for(int col = 0; col < states; col++){
                         (*A)[row][col] = numbers[row*states+col];
                     }
                 }

                 break;
             }
             case 1://Resize emission matrix and fill in data
             {
                 measurments = numbers[1];
                 numbers.erase(numbers.begin());
                 numbers.erase(numbers.begin());

                 (*B).resize(states);
                 for(int i = 0; i < states; i++){
                     (*B)[i].resize(measurments);
                 }

                 for(int row = 0; row < states; row++){
                     for(int col = 0; col < measurments; col++){
                         (*B)[row][col] = numbers[row*measurments+col];
                     }
                 }

                 break;
             }
             case 2://Resize initialization vector and fill in data
             {
                 (*pi).resize(states);
                 numbers.erase(numbers.begin());
                 numbers.erase(numbers.begin());

                 for(int row = 0; row < states; row++){
                     (*pi)[row] = numbers[row];
                 }
                 break;
             }
             case 3://Resize sequence vector and fill in data
             {
                 tMax = numbers[0];
                 numbers.erase(numbers.begin());
                 (*seq).resize(tMax);

                 for(int row = 0; row < tMax; row++){
                     (*seq)[row] = numbers[row];
                 }
                 break;
             }
         }
         matrix++;
     }
}

double alphaLDPass(std::vector<std::vector<double>>* A,std::vector<std::vector<double>>* B, std::vector<double>* pi, std::vector<double>* seq, std::vector<std::vector<double>>* alpha){
    //Calculates the probability of emitting an observed sequence(seq) at a given time step(t=tMax-1)

    //Fill the first row of alpha
    //alpha_i(1) = pi_i*b_i(O_1)
    for(int i = 0; i < states; i++){
        (*alpha)[0][i] = elnproduct(eln((*pi)[i]),eln((*B)[i][(*seq)[0]]));
    }

    //Fill alpha using dynamic programming
    //Alpha_i(t) = sum{j:0->N}[alpha_t-1(j)*a_ji]*b_i(O_t)
    for (int t = 1; t < tMax; t++) {
        for(int i = 0; i < states; i++) {
            double sum = LOGZERO;
            for (int j = 0; j < states; j++) {
                double lnprod = elnproduct((*alpha)[t-1][j],eln((*A)[j][i]));
                sum = elnsum(sum,lnprod);
            }

            double secondFactor = eln((*B)[i][((*seq)[t])]);
            (*alpha)[t][i] = elnproduct(sum,secondFactor);
        }
    }

    //Extract solution
    //P(O_1:tMax) = sum{i:0->N}[alpha[tMax-1][i]]
    double ans = LOGZERO;
    for(int i = 0; i < states; i++){
        ans =  elnsum(ans,(*alpha)[tMax-1][i]);
    }

    return ans;
}

void betaLDPass(std::vector<std::vector<double>>* A, std::vector<std::vector<double>>* B, std::vector<double>* seq, std::vector<std::vector<double>>* beta){
    //Calculates the probability of emitting an observed sequence(seq) at a given time step(t=tMax-1)

    //Fill the first row of beta
    //beta_i(1) = pi_i*b_i(O_1)
    for(int i = 0; i < states; i++){
        (*beta)[tMax-1][i] = 0;
    }

    //Fill beta using dynamic programming
    //Beta_i(t) = sum{j:0->N}[beta_t-1(j)*a_ji]*b_i(O_t)
    for (int t = tMax-2; t >= 0; t--) {
        for(int i = 0; i < states; i++) {
            (*beta)[t][i] = LOGZERO;
            for (int j = 0; j < states; j++) {
                double secondFactor = elnproduct(eln((*B)[j][(*seq)[t+1]]),(*beta)[t+1][j]);
                double finalProduct = elnproduct(eln((*A)[i][j]),secondFactor);
                (*beta)[t][i] = elnsum((*beta)[t][i],finalProduct);
            }
        }
    }
}

double getGammaLDDenominator(std::vector<std::vector<double>>* alpha){
    //Calculates the denominator for the gamma and digamma calculations
    double denom = LOGZERO;
    for(int i = 0; i < states; i++){
        denom = elnsum(denom,(*alpha)[tMax-1][i]);
    }
    return denom;
}

void gammaLDFun(std::vector<std::vector<double>>* alpha, std::vector<std::vector<double>>* beta, std::vector<std::vector<double>>* gamma){
    //Calculates gamma
    double denom = getGammaLDDenominator(alpha);

    for(int t = 0; t < tMax; t++){
        for(int i = 0; i < states; i++){
            (*gamma)[t][i] = elnproduct(elnproduct(((*alpha)[t][i]),((*beta)[t][i])),-denom);
        }
    }
}

void diGammaLDFun(std::vector<std::vector<double>>* A, std::vector<std::vector<double>>* B, std::vector<double>* seq, std::vector<std::vector<double>>* alpha, std::vector<std::vector<double>>* beta, std::vector<std::vector<std::vector<double>>>* diGamma){
    //Calculates digamma
    double denom = getGammaLDDenominator(alpha);

    for(int t = tMax-2; t >= 0; t--){
        for(int i = 0; i < states; i++){
            for(int j = 0; j < states; j++){
                double alphaVal = (*alpha)[t][i];
                double aVal = eln((*A)[i][j]);
                double bVal = eln((*B)[j][((*seq)[t+1])]);
                double betaVal = (*beta)[t+1][j];
                double numerator = elnproduct(elnproduct(alphaVal,aVal),elnproduct(bVal,betaVal));
                (*diGamma)[t][i][j] = elnproduct(numerator,-denom);
            }
        }
    }
}

void reEstimateParamsLD(std::vector<std::vector<double>>* A, std::vector<std::vector<double>>* B, std::vector<double>* pi, std::vector<double>* seq, std::vector<std::vector<double>>* gamma, std::vector<std::vector<std::vector<double>>>* diGamma){
    //Reestimates the model parameters
    //Transition matrix
    for(int i = 0; i < states; i++){
        double denominator = LOGZERO;
        for(int t = 0; t <= tMax-2; t++){
            denominator = elnsum(denominator,(*gamma)[t][i]);
        }
        for(int j = 0; j < states; j++){
            double numerator = LOGZERO;
            for(int t = 0; t <= tMax-2; t++){
                numerator = elnsum(numerator,(*diGamma)[t][i][j]);
            }
            (*A)[i][j] = eexp(elnproduct(numerator,-denominator));//Comment out to skip transitionmatrix
        }
    }


    //Emission matrix
    for(int i = 0; i < states; i++){
        for(int k = 0; k < measurments; k++){
            double denominator = LOGZERO;
            double numerator = LOGZERO;
            for (int t = 0;t <= tMax-2; t++){
                denominator = elnsum(denominator,(*gamma)[t][i]);
                if((*seq)[t]==k){
                    numerator = elnsum(numerator,(*gamma)[t][i]);
                }
            }
            if (denominator == LOGZERO){
                std::cerr << "WARNING: DIVISION BY ZERO";
                continue;
            }
            (*B)[i][k] = eexp(elnproduct(numerator,-denominator));
        }
    }

    //Initialization vector
    for(int i = 0; i < states; i++){
        (*pi)[i] = eexp((*gamma)[0][i]);
    }

}

void outputMatrix(std::vector<std::vector<double>> matrix){
    //Outputs a Matrix's dimensions followed by its values
    std::cout.precision(20);
    std::cout << matrix.size() << " " << matrix[0].size() << " ";
    for (int r = 0; r < (int) matrix.size(); r++) {
        for (int c = 0; c < (int) matrix[0].size(); c++) {
            std::cout << round(1000000000000000000*matrix[r][c])/1000000000000000000 << + " ";//Round very small values for convenience
        }
    }
    std::cout << std::endl;
}

//-----Log Helper Functions----
//Scaling inspired by:
//http://bozeman.genome.washington.edu/compbio/mbt599_2006/hmm_scaling_revised.pdf
double eexp(double x){
    //Calculates e^x
    if(x == LOGZERO){
        return 0;
    }
    return exp(x);
}

double eln(double x){
    //Calculates ln(x)
    if(x == 0){
        return LOGZERO;
    }
    return log(x);
}

double elnsum(double elnx, double elny){
    //Calculates ln(x)+ln(Y)
    if (elnx == LOGZERO){
        return elny;
    }
    else if(elny == LOGZERO){
        return elnx;
    }
    else if(elnx > elny){
        return elnx + eln(1+exp(elny-elnx));
    }
    else{
        return elny + eln(1+exp(elnx-elny));
    }
}

double elnproduct(double elnx,double elny){
    //Calculates ln(x)*ln(Y)
    if (elnx == LOGZERO or elny == LOGZERO){
        return LOGZERO;
    }
    else{
        return elnx+elny;
    }
}

//-----Debbugging Functions----
void displayMatrix(std::vector<std::vector<double>> matrix){
    //Displays a specified matrix in stderr
    for(int r = 0; r < (int)matrix.size(); r++){
        for(int c = 0; c < (int)matrix[0].size(); c++){
            std::cerr << matrix[r][c] << " ";
        }
        std::cerr << std::endl;
    }
    std::cerr << std::endl;
}

void displayMatrixLD(std::vector<std::vector<double>> matrix){
    //Displays a specified log domain matrix in stderr
    for(int r = 0; r < (int)matrix.size(); r++){
        for(int c = 0; c < (int)matrix[0].size(); c++){
            std::cerr << eexp(matrix[r][c]) << " ";
        }
        std::cerr << std::endl;
    }
    std::cerr << std::endl;
}

void displayVector(std::vector<double> vector){
    //Displays a specified vector in stderr
    for(int r = 0; r < (int)vector.size(); r++){
        std::cerr << vector[r] << " ";
    }
    std::cerr << std::endl << std::endl;
}

//--------Main Function--------
int main() {
    std::vector<std::vector<double>> A;//Declare transition matrix a_ij
    std::vector<std::vector<double>> B;//Declare emission matrix b_ik
    std::vector<double> pi;//Declare initialization vector
    std::vector<double> seq;//Declare sequence vector

    readInput(&A,&B,&pi,&seq);

    //Create the alpha, beta, gamma and di-gamma matrices
    std::vector<std::vector<double>> alphaLD(tMax);
    std::vector<std::vector<double>> betaLD(tMax);
    std::vector<std::vector<double>> gammaLD(tMax);
    std::vector<std::vector<std::vector<double>>> diGammaLD(tMax);
    for(int t = 0; t < tMax; t++){
        alphaLD[t].resize(states);
        betaLD[t].resize(states);
        gammaLD[t].resize(states);
        diGammaLD[t].resize(states);

        for(int i = 0; i < states; i++){
            diGammaLD[t][i].resize(states);
        }
    }

    double prob = 0;
    double newProb = -1;
    while (std::abs(prob-newProb)>0.01){//Baum Welch
        prob = newProb;
        alphaLDPass(&A, &B, &pi, &seq, &alphaLD);
        betaLDPass(&A, &B, &seq, &betaLD);
        gammaLDFun(&alphaLD, &betaLD, &gammaLD);
        diGammaLDFun(&A, &B, &seq, &alphaLD, &betaLD, &diGammaLD);

        reEstimateParamsLD(&A, &B, &pi, &seq, &gammaLD, &diGammaLD);

        newProb = alphaLDPass(&A, &B, &pi, &seq, &alphaLD);
    }

    //Output the answer
    outputMatrix(A);
    outputMatrix(B);

    return 0;
}
