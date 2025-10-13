/*
   This file is part of the RELXILL model code.

   RELXILL is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   any later version.

   RELXILL is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.
   For a copy of the GNU General Public License see
   <http://www.gnu.org/licenses/>.

    Copyright 2022 Thomas Dauser, Remeis Observatory & ECAP
*/
#ifndef RELXILL_ERRORHANDLING_H_
#define RELXILL_ERRORHANDLING_H_

#include <exception>
#include <string>


/**
 * @brief Exception thrown when a status variable indicates EXIT_FAILURE
 *
 * This exception is thrown when a previous function call exits with failure status.
 * It captures the function name where the error occurred and provides a descriptive
 * error message following C++ standard exception handling patterns.
 */
class RelxillStatusException : public std::exception
{
public:
    RelxillStatusException() = default;

    /**
     * @brief Construct exception with function name
     * @param func_name Name of the function where the error occurred
     */
    explicit RelxillStatusException(const char* func_name)
    {
        m_msg += "error in function ";
        m_msg += func_name;
        m_msg += " from previous call exiting with failure status";
    }

    /**
     * @brief Get the error message
     * @return C-string containing the error message
     */
    [[nodiscard]] const char* what() const noexcept override
    {
        return m_msg.c_str();
    }

private:
    std::string m_msg{"*** relxill-error: "};
};

/**
 * @brief Check status and throw exception if it indicates failure
 *
 * This function checks if the status variable equals EXIT_FAILURE and throws
 * a RelxillStatusException with the calling function's name.
 *
 * @param func_name Name of the calling function (typically __func__)
 * @param status Pointer to status variable to check
 * @throws RelxillStatusException if *status == EXIT_FAILURE
 */
void check_status_and_throw(const char* func_name, int status);

// Macro to automatically capture the function name
#define CHECK_STATUS_AND_THROW(status) check_status_and_throw(__func__, status);

#endif /* RELXILL_ERRORHANDLING_H_ */