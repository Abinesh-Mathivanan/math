#include <vector>
#include <random>
#include <cmath>
#include <memory>

class GroupRPO {
private:
    int num_agents;
    int state_dim;
    int action_dim;
    double learning_rate;
    double clip_ratio;
    std::vector<std::vector<double>> policies;
    
public:
    GroupRPO(int num_agents, int state_dim, int action_dim, 
             double learning_rate = 0.01, double clip_ratio = 0.2) 
        : num_agents(num_agents)
        , state_dim(state_dim)
        , action_dim(action_dim)
        , learning_rate(learning_rate)
        , clip_ratio(clip_ratio) {
        policies.resize(num_agents);
        for (auto& policy : policies) {
            policy.resize(state_dim * action_dim, 0.0);
        }
    }

    std::vector<double> computePolicy(int agent_id, const std::vector<double>& state) {
        std::vector<double> action_probs(action_dim);
        double sum = 0.0;
        
        for (int i = 0; i < action_dim; i++) {
            double logit = 0.0;
            for (int j = 0; j < state_dim; j++) {
                logit += policies[agent_id][i * state_dim + j] * state[j];
            }
            action_probs[i] = std::exp(logit);
            sum += action_probs[i];
        }

        for (double& prob : action_probs) {
            prob /= sum;
        }
        
        return action_probs;
    }

    void update(const std::vector<std::vector<double>>& states,
                const std::vector<std::vector<int>>& actions,
                const std::vector<std::vector<double>>& advantages) {
        
        for (int agent = 0; agent < num_agents; agent++) {
            std::vector<double> policy_grad(state_dim * action_dim, 0.0);
            
            for (size_t t = 0; t < states.size(); t++) {
                auto action_probs = computePolicy(agent, states[t]);
                double ratio = action_probs[actions[t][agent]];
                
                double clipped_ratio = std::clamp(ratio,
                    1.0 - clip_ratio,
                    1.0 + clip_ratio);
                
                double surrogate = std::min(
                    ratio * advantages[t][agent],
                    clipped_ratio * advantages[t][agent]
                );
                
                for (int i = 0; i < state_dim * action_dim; i++) {
                    policies[agent][i] += learning_rate * surrogate * 
                        (states[t][i % state_dim]);
                }
            }
        }
    }

    std::vector<int> sampleActions(const std::vector<double>& state) {
        std::vector<int> actions(num_agents);
        std::random_device rd;
        std::mt19937 gen(rd());
        
        for (int agent = 0; agent < num_agents; agent++) {
            auto probs = computePolicy(agent, state);
            std::discrete_distribution<> dist(probs.begin(), probs.end());
            actions[agent] = dist(gen);
        }
        
        return actions;
    }
};