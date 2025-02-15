#include <vector>
#include <random>
#include <cmath>

class CNN {
private:
    struct Layer {
        std::vector<std::vector<std::vector<float>>> filters;
        std::vector<std::vector<std::vector<float>>> output;
        std::vector<float> bias;
    };

    std::vector<Layer> layers;
    float learning_rate;

    float relu(float x) {
        return std::max(0.0f, x);
    }

    float relu_derivative(float x) {
        return x > 0 ? 1.0f : 0.0f;
    }

    std::vector<std::vector<float>> convolve(const std::vector<std::vector<float>>& input,
                                           const std::vector<std::vector<float>>& kernel) {
        int output_height = input.size() - kernel.size() + 1;
        int output_width = input[0].size() - kernel[0].size() + 1;
        std::vector<std::vector<float>> output(output_height, std::vector<float>(output_width));

        for (int i = 0; i < output_height; ++i) {
            for (int j = 0; j < output_width; ++j) {
                float sum = 0;
                for (int k = 0; k < kernel.size(); ++k) {
                    for (int l = 0; l < kernel[0].size(); ++l) {
                        sum += input[i + k][j + l] * kernel[k][l];
                    }
                }
                output[i][j] = sum;
            }
        }
        return output;
    }

public:
    CNN(float lr = 0.01) : learning_rate(lr) {}

    void add_layer(int num_filters, int filter_size) {
        Layer layer;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<float> dist(0.0f, 0.1f);

        for (int i = 0; i < num_filters; ++i) {
            std::vector<std::vector<float>> filter(filter_size, std::vector<float>(filter_size));
            for (int j = 0; j < filter_size; ++j) {
                for (int k = 0; k < filter_size; ++k) {
                    filter[j][k] = dist(gen);
                }
            }
            layer.filters.push_back(filter);
            layer.bias.push_back(0.0f);
        }
        layers.push_back(layer);
    }

    std::vector<std::vector<float>> forward(const std::vector<std::vector<float>>& input) {
        std::vector<std::vector<float>> current_input = input;

        for (auto& layer : layers) {
            std::vector<std::vector<std::vector<float>>> layer_output;
            for (size_t i = 0; i < layer.filters.size(); ++i) {
                auto conv_output = convolve(current_input, layer.filters[i]);
                for (auto& row : conv_output) {
                    for (auto& val : row) {
                        val = relu(val + layer.bias[i]);
                    }
                }
                layer_output.push_back(conv_output);
            }
            layer.output = layer_output;
            if (layer_output.size() == 1) {
                current_input = layer_output[0];
            }
        }
        return current_input;
    }

    void backward(const std::vector<std::vector<float>>& input,
                 const std::vector<std::vector<float>>& target) {
        std::vector<std::vector<float>> error;
        for (int i = 0; i < target.size(); ++i) {
            std::vector<float> row;
            for (int j = 0; j < target[0].size(); ++j) {
                row.push_back(target[i][j] - layers.back().output[0][i][j]);
            }
            error.push_back(row);
        }

        for (int l = layers.size() - 1; l >= 0; --l) {
            for (size_t f = 0; f < layers[l].filters.size(); ++f) {
                std::vector<std::vector<float>> filter_gradient(layers[l].filters[f].size(),
                    std::vector<float>(layers[l].filters[f][0].size(), 0.0f));

                for (int i = 0; i < error.size(); ++i) {
                    for (int j = 0; j < error[0].size(); ++j) {
                        float gradient = error[i][j] * relu_derivative(layers[l].output[f][i][j]);
                        layers[l].bias[f] += learning_rate * gradient;

                        for (int k = 0; k < layers[l].filters[f].size(); ++k) {
                            for (int m = 0; m < layers[l].filters[f][0].size(); ++m) {
                                filter_gradient[k][m] += gradient * input[i + k][j + m];
                            }
                        }
                    }
                }

                for (size_t i = 0; i < layers[l].filters[f].size(); ++i) {
                    for (size_t j = 0; j < layers[l].filters[f][0].size(); ++j) {
                        layers[l].filters[f][i][j] += learning_rate * filter_gradient[i][j];
                    }
                }
            }
        }
    }
};