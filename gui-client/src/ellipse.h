#include <string>
#include <vector>

void calculateEllipsesParameters(std::vector<std::string>& out, std::vector<std::vector<float>>& x1El,
                                 std::vector<std::vector<float>>& x2El, std::vector<std::vector<float>>& yEl1,
                                 std::vector<std::vector<float>>& yEl2, const std::vector<float>& x,
                                 const std::vector<float>& y, const std::vector<float>& targets);
void plotEllipse(float ra, float rb, float alpha, std::vector<float>& x1, std::vector<float>& x2,
                 std::vector<float>& y1, std::vector<float>& y2, float mx, float my);
float inEllipse(float ra, float rb, float alpha, const std::vector<float>& x, const std::vector<float>& y, float mx,
                float my);
float sum(const std::vector<float>& v);
float disp(const std::vector<float>& v);
float K(const std::vector<float>& vx, const std::vector<float>& vy);
std::pair<std::vector<float>, std::vector<float>> calc_alpha_beta(const std::vector<std::vector<float>>& x, const std::vector<std::vector<float>>& y, const float target);
std::pair<std::vector<float>, std::vector<float>> calc_cm(const std::vector<std::vector<float>>& x,
	const std::vector<std::vector<float>>& y);