#include "input_models.hpp"
#include <algorithm>
#include <fstream>
#include <regex>
#include <stdexcept>

namespace gasp2::details {

// Quickly inspect the XML input to determine which surface and reaction models
// are requested and whether any reaction definitions are present.
[[nodiscard]] InputModels
read_input_models(const std::filesystem::path &input_filename) {
  if (!std::filesystem::exists(input_filename)) {
    throw std::runtime_error("Input file does not exist: " +
                             input_filename.string());
  }

  std::ifstream in(input_filename);
  // Read entire file into string using iterator-based constructor
  // First iterator reads from file, second is end-of-stream sentinel
  std::string xml((std::istreambuf_iterator<char>(in)),
                  std::istreambuf_iterator<char>());

  InputModels info;

  // Check if XML contains any <reaction> tags (case-insensitive, allows
  // whitespace)
  std::regex reaction_tag_re{R"(<\s*reaction\b)", std::regex::icase};
  info.has_reactions = std::regex_search(xml, reaction_tag_re);

  // Container to store regex match results (capture groups)
  std::smatch fm;
  // Extract surface_model attribute from <gasp2> tag
  // Pattern: <gasp2 [any attrs] surface_model="(value)" [any attrs]>
  std::regex surface_model_re{R"(<gasp2\s+[^>]*surface_model=\"([^\"]+)\")",
                              std::regex::icase};
  if (std::regex_search(xml, fm, surface_model_re)) {
    // Extract attribute value and convert to lowercase for comparison
    std::string val = fm[1].str();
    std::transform(val.begin(), val.end(), val.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    // Map string value to corresponding enum
    if (val == "catalysis")
      info.surface_model = SurfaceModel::Catalysis;
    else if (val == "non_catalytic")
      info.surface_model = SurfaceModel::NonCatalytic;
    else if (val == "ablation")
      info.surface_model = SurfaceModel::Ablation;
    else
      throw std::runtime_error("Unknown surface model: " + val);
  } else {
    throw std::runtime_error("surface_model attribute missing");
  }

  std::regex model_re{R"(<reactions\s+[^>]*model=\"([^\"]+)\")",
                      std::regex::icase};
  if (std::regex_search(xml, fm, model_re)) {
    std::string val = fm[1].str();
    std::transform(val.begin(), val.end(), val.begin(),
                   [](unsigned char c) { return std::tolower(c); });
    if (val == "gamma_model")
      info.model = ReactionsModel::GammaModel;
    else if (val == "finite_rate" || val == "finite_rates")
      info.model = ReactionsModel::FiniteRates;
    else
      throw std::runtime_error("Unknown reactions model: " + val);
  }

  return info;
}

} // namespace gasp2::details
