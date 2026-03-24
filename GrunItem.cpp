#include "GrunItem.h"

const GrunItem GrunItem::INVALID("INVALID_ITEM_SENTINEL", -1);

GrunItem::GrunItem(	std::string name,
					int			id,
					std::string relationship,
					std::string relComment,
					std::string quantityFormula, 
					std::string units,
					std::string primaryLabourFormula
				)
				:	_itemName(std::move(name)),
					_libraryId(id),
					_itemQuantityFormula(std::move(quantityFormula)),
					_itemQuantityUnits(std::move(units)),
					_itemPrimaryLabourFormula(std::move(primaryLabourFormula))
{
	// on GrunItem creation, create a CoreValues object with default values and push it into the _itemCoreValues vector
	addRelationshipValues(relationship, relComment, 0.0, 0.0, SpatialExponentValue::None);
}

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
 * @brief Adds a new set of RelationshipValues to the GrunItem. All arguments are optional.
 * @param relationship				(std::string)			- the GrunItem's relationship string (default: empty)
 * @param relComment				(std::string)			- this relationship's optional comment string (default: empty)
 * @param itemQuantity				(double)				- gets calculated later (default: 0.0)
 * @param spatialQuantity			(double)				- gets calculated later (default: 0.0)
 * @param spatialUnit				(SpatialExponentValue)	- gets calculated later (default: SpatialExponentValue::None)
 * @return	(int) total number of relationships the GrunItem has after this new addition
 */
size_t GrunItem::addRelationshipValues(	const std::string& relationship,
										const std::string& relComment,
										const double& itemQuantity,
										const double& spatialQuantity,
										const SpatialExponentValue& spatialUnit
									)
{
	// this method is *strictly* for ADDING a new RelationshipValues set to the GrunItem
	RelationshipValues	newValues;
	newValues.relationship				= relationship;
	newValues.relComment				= relComment;
	newValues.itemQuantity				= itemQuantity;
	newValues.spatialQuantity			= spatialQuantity;
	newValues.spatialUnit				= spatialUnit;

	this->_itemCoreValues.push_back(newValues);
	return this->_itemCoreValues.size();
}

/** 
 * @brief Removes a RelationshipValues from the _itemCoreValues vector by index value
 */
size_t GrunItem::rmRelationshipValues(const size_t index)
{
	// zero-check
	if (index >= this->_itemCoreValues.size())
	{
		return 0;
	}
	if (index >= 0 && index < this->_itemCoreValues.size())
	{
		this->_itemCoreValues.erase(this->_itemCoreValues.begin() + index);
		return this->_itemCoreValues.size();
	}
	return size_t();
}

/**
 * @brief UPDATES an existing relationship in an existing GrunItem
 * 
 * @param index 		(size_t)		the location of the relationship in the GrunItem's coreValue vector
 * @param relationship  (std::string)	the new relationship string
 * @param relComment 	(std::string)	the new relComment string
 * @return true  if the relationship is updated successfully
 * @return false if the provided index is out of bounds
 */
GrunItem& GrunItem::updateRelationshipValue(const size_t index, const std::string relationship, const std::string relComment)
{
	// zero-check
	if (index >= this->_itemCoreValues.size())
	{
		return const_cast<GrunItem&>(GrunItem::INVALID);
	}
	if (index >= 0 && index < this->_itemCoreValues.size())
	{
		this->_itemCoreValues[index].relationship 	= relationship;
		this->_itemCoreValues[index].relComment		= relComment;
	}

	if (onValueChange) { onValueChange(*this); }
	return *this;
}

/**
 * @brief Converts a GrunItem's time-typed member to user-friendly date/time string, returned in optional format
 * @param member 	- the time-typed member in the GrunItem we want to retrieve (required)
 * @param format	- std::string that defines the format the time should be shown in, look up std::strftime() at https://en.cppreference.com/w/cpp/chrono/c/strftime.html to see the tokens you can use (default: "%Y%m%d %H:%M:%S")
 * @return std::string	- The formatted datetime string, NULL" if the item's attributes have never been calculated, or error msg if an error is encountered
 */
std::string	GrunItem::getCalculatedTimeString(const std::chrono::system_clock::time_point& member, const std::string& format) 
{ 
	/* 
		dev-note: this function uses some archaic looking C-style shit, because apparently if we want to use the "format" argument to 
		allow us to customize the way the timestamp is displayed in output, C++20 doesnt have a way to do this, so instead we have to 
		convert everything back into ancient C types and use buffers and shit to make it compile.
	*/

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
