
# Custom logger type that includes timestamps
struct TimeStampedLogger <: AbstractLogger
    io::IO  # File handle to write logs to
    min_level::Logging.LogLevel  # Minimum log level to record
end

# Define the behavior for handling log messages
function Logging.handle_message(logger::TimeStampedLogger, level::Logging.LogLevel, message, _module, group, id, file, line; kwargs...)
    if level >= logger.min_level  # Only log messages that meet the minimum log level
        timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
        println(logger.io, "[$timestamp] ", Logging.level_string(level), " ", message)
        flush(logger.io)  # Ensure the message is written immediately
    end
end
# Define the minimum log level method
function Logging.min_enabled_level(logger::TimeStampedLogger)
    return logger.min_level
end


# Function to configure timestamped file logger
function configure_file_logger(logfile::String)
    open(logfile, "a") do io  # Open the file in append mode
        global_logger(TimeStampedLogger(io, Logging.Error))  # Set logger to record errors and above
    end
end

# Configure the logger to write to "error_log.txt"

# Example of a try-catch block where an error is logged
