#include "GrunItem.h"

/**
 * @brief Returns a SpatialExponentValue as a string for friendly output
 * @param exponent	- the SpatialExponentValue to return
 * @return std::string	- The value of the SpatialExponentValue as a std::string, else "UNKNOWN" if something goes wrong
 */
std::string spatialExponentValueToString(SpatialExponentValue exponent)
{
	switch (exponent)
    {
        case SpatialExponentValue::None:	return "None";
        case SpatialExponentValue::Linear:		return "Linear";
        case SpatialExponentValue::Area:		return "Area";
        case SpatialExponentValue::Volume:		return "Volume";
        default:					            return "UNKNOWN";
    }
}

/**
 * @brief Converts a GrunItem's time-typed member to user-friendly date/time string, returned in optional format
 * @param member 	- the time-typed member in the GrunItem we want to retrieve (required)
 * @param format	- std::string that defines the format the time should be shown in (default: "%Y%m%d %H:%M:%S")
 * @return std::string	- The formatted datetime string, NULL" if the item's attributes have never been calculated, or error msg if an error is encountered
 */
std::string	GrunItem::getCalculatedTimeString(const std::chrono::system_clock::time_point& member, const std::string& format) 
{ 
	// dev-note: this function uses some archaic looking C-style shit, because apparently if we want to use the "format" argument to 
	// allow us to customize the way the timestamp is displayed in output, C++20 doesnt have a way to do this, so instead we have to 
	// convert everything back into ancient C types and use buffers and shit to make it compile.

	// zero-check: if the p_member's value is in the default state, it means the GrunItem's attributes have never been calculated, so return the LKGWCalculatedTime string as "NULL"
	if (std::chrono::system_clock::to_time_t(member) == 0) { return std::string("NULL"); }
	try
	{
		// convert time_point to std::time_t
		const std::time_t t_c = std::chrono::system_clock::to_time_t(member);
		// convert time_t to local time tm structure
		std::tm* tm_local = std::localtime(&t_c);
		if (!tm_local) { return std::string("Time conversion error"); }

		std::string buffer(128,'\0');
		size_t size = buffer.size();
		size_t written = 0;	
		while (true)
		{
			written = std::strftime(buffer.data(), size, format.c_str(), tm_local);
			if (written > 0)
			{
				buffer.resize(written);
				return buffer;
			}
			else if (written == 0 && size == 0)
			{
				if (buffer.size() > 1024)
				{
					return std::string("Time formatting error: buffer limit exceeded");
				}
				size *= 2;
				buffer.resize(size);
			}
			else
			{
				return std::string("Time formmating failed.");
			}
		}
	}
	catch(const std::exception& e)
	{
		// catch exceptions that could arise from bad format strings
		return std::format("Caught exception, the datetime format string is invalid: '{}'", e.what());
	}
}