#include <vector>
#include <random>
#include <algorithm>

class HiddenMarkovModel {
private:
    std::vector<std::vector<double>> transitionMatrix;
    std::vector<std::vector<double>> emissionMatrix;
    std::vector<double> initialState;
    int numStates;
    int numObservations;

public:
    HiddenMarkovModel(int states, int observations) : 
        numStates(states), 
        numObservations(observations) {
        transitionMatrix.resize(states, std::vector<double>(states));
        emissionMatrix.resize(states, std::vector<double>(observations));
        initialState.resize(states);
    }

    void setTransitionMatrix(const std::vector<std::vector<double>>& matrix) {
        transitionMatrix = matrix;
    }

    void setEmissionMatrix(const std::vector<std::vector<double>>& matrix) {
        emissionMatrix = matrix;
    }

    void setInitialState(const std::vector<double>& initial) {
        initialState = initial;
    }

    std::vector<int> viterbi(const std::vector<int>& observations) {
        std::vector<std::vector<double>> dp(observations.size(), std::vector<double>(numStates));
        std::vector<std::vector<int>> backpointer(observations.size(), std::vector<int>(numStates));

        for (int i = 0; i < numStates; i++) {
            dp[0][i] = initialState[i] * emissionMatrix[i][observations[0]];
        }

        for (int t = 1; t < observations.size(); t++) {
            for (int j = 0; j < numStates; j++) {
                double maxProb = 0.0;
                int maxState = 0;

                for (int i = 0; i < numStates; i++) {
                    double prob = dp[t-1][i] * transitionMatrix[i][j] * emissionMatrix[j][observations[t]];
                    if (prob > maxProb) {
                        maxProb = prob;
                        maxState = i;
                    }
                }

                dp[t][j] = maxProb;
                backpointer[t][j] = maxState;
            }
        }

        std::vector<int> path(observations.size());
        double maxProb = 0.0;
        int maxState = 0;

        for (int i = 0; i < numStates; i++) {
            if (dp[observations.size()-1][i] > maxProb) {
                maxProb = dp[observations.size()-1][i];
                maxState = i;
            }
        }

        path[observations.size()-1] = maxState;
        for (int t = observations.size()-2; t >= 0; t--) {
            path[t] = backpointer[t+1][path[t+1]];
        }

        return path;
    }

    std::vector<double> forward(const std::vector<int>& observations) {
        std::vector<std::vector<double>> alpha(observations.size(), std::vector<double>(numStates));

        for (int i = 0; i < numStates; i++) {
            alpha[0][i] = initialState[i] * emissionMatrix[i][observations[0]];
        }

        for (int t = 1; t < observations.size(); t++) {
            for (int j = 0; j < numStates; j++) {
                double sum = 0.0;
                for (int i = 0; i < numStates; i++) {
                    sum += alpha[t-1][i] * transitionMatrix[i][j];
                }
                alpha[t][j] = sum * emissionMatrix[j][observations[t]];
            }
        }

        std::vector<double> result(numStates);
        for (int i = 0; i < numStates; i++) {
            result[i] = alpha[observations.size()-1][i];
        }
        return result;
    }
};