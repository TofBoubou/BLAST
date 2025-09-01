#pragma once
#include "gasp2/types.hpp"
#include <filesystem>

namespace gasp2::details {

/// Minimal information about the selected models extracted from the input.
struct InputModels {
  ReactionsModel model{ReactionsModel::GammaModel};
  SurfaceModel surface_model{SurfaceModel::Catalysis};
  bool has_reactions{false}; ///< Whether any <reaction> tags are present.
};

/// Quickly scan the XML file to determine which models are requested and
/// whether any reaction blocks are present.
[[nodiscard]] InputModels
read_input_models(const std::filesystem::path &input_filename);

} // namespace gasp2::details
