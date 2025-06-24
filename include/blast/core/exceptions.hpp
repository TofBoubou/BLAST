#pragma once
#include <stdexcept>
#include <string>
#include <string_view>
#include <source_location>
#include <format>

namespace blast::core {

class BlastException : public std::exception {
private:
    std::string message_;
    std::source_location location_;
    
public:
    explicit BlastException(std::string message, std::source_location location = std::source_location::current())
        : message_(message), location_(location) {}
    
    template<typename... Args>
    explicit BlastException(std::string_view format_str, std::source_location location, Args&&... args) 
        noexcept : location_(location) {
        try {
            message_ = std::format(format_str, std::forward<Args>(args)...);
        } catch (...) {
            message_ = "Exception formatting failed";
        }
    }
    
    [[nodiscard]] const char* what() const noexcept override {return message_.c_str();}
    
    [[nodiscard]] auto location() const noexcept -> const std::source_location& {return location_;}
    
    [[nodiscard]] auto message() const noexcept -> const std::string& {return message_;}

    [[nodiscard]] auto full_message() const -> std::string {
        return std::format(
            "{} [{}:{}:{}]",
            message_,
            location_.file_name(),
            location_.line(),
            location_.function_name()
        );
    }
};


class ConfigurationError : public BlastException {
public:
    explicit ConfigurationError(std::string_view message, std::source_location location = std::source_location::current()) 
        : BlastException(std::format("Configuration Error: {}", message), location) {}
    
    template<typename... Args>
    explicit ConfigurationError(std::string_view format_str, std::source_location location, Args&&... args) 
        : BlastException(std::format("Configuration Error: {}", std::format(format_str, std::forward<Args>(args)...)), location) {}
};

class FileError : public BlastException {
private:
    std::string filename_;
    
public:
    explicit FileError(std::string_view message, std::string filename, std::source_location location = std::source_location::current()) 
        : BlastException(std::format("File Error ({}): {}", filename, message), location), filename_(std::move(filename)) {}
    
    [[nodiscard]] auto filename() const noexcept -> const std::string& { return filename_; }
};

class ValidationError : public ConfigurationError {
private:
    std::string field_name_;
    
public:
    explicit ValidationError(std::string_view field_name, std::string_view message, std::source_location location = std::source_location::current()) 
        : ConfigurationError(std::format("Field '{}': {}", field_name, message), location), field_name_(field_name) {}
    
    [[nodiscard]] auto field_name() const noexcept -> const std::string& { return field_name_; }
};

class GridError : public core::BlastException {
public:
    explicit GridError(std::string_view message, 
                      std::source_location location = std::source_location::current()) noexcept
        : BlastException(std::format("Grid Error: {}", message), location) {}
};

class TransformError : public core::BlastException {
public:
    explicit TransformError(std::string_view message,
                           std::source_location location = std::source_location::current()) noexcept
        : BlastException(std::format("Transform Error: {}", message), location) {}
};

}