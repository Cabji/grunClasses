#include <algorithm>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <print>
#include <ranges>
#include <regex>
#include <source_location>
#include <stdexcept>
#include <sstream>
#include "cabjiHelpers.h"
#include "GrunObject.h"

const double PI = std::acos(-1.0);

const	std::string	GO_LINEAL_TOKENS					= "LWDCR";
const	std::string	GO_AREA_TOKENS						= "A";
const	std::string	GO_VOLUME_TOKENS					= "V";
const	std::string	GO_SPATIAL_SIGNIFICANT_OPERATORS	= "*";

// regex pattern strings that are stored as constants
// this is mostly to store these regexes in 1 convenient place in in the source code so you don't have to go searching for them.
const	std::regex	REGEX_GI_BASEEXPR_SIG_TOKENS_AND_OPS(R"(([^LWDAVCR\/*]))");;
const	std::regex	REGEX_GO_ALL_TOKENS(R"([LWDCRAV])");
const	std::regex	REGEX_GO_LINEAL_TOKENS(R"([LWDCR]+)");
const	std::regex	REGEX_GO_AREA_TOKENS(R"([A]+)");
const	std::regex	REGEX_GO_VOLUME_TOKENS(R"([V]+)");
const	std::regex	REGEX_SHN_TO_PEDMAS_0_WRAP_ALL_IN_PARENTHESES(R"((.+))");
const	std::regex	REGEX_SHN_TO_PEDMAS_1_EXPLICIT_COMBINE_OPERATOR(R"(([+\-]))");
const	std::regex	REGEX_SHN_TO_PEDMAS_2_NUM_FACTOR_AND_GO_TOKEN(R"(([\d.]*)([LWDAVCR]))");
const	std::regex	REGEX_SHN_TO_PEDMAS_3_AT_OPERATOR(R"((@)([\d.]+))");
const	std::regex	REGEX_SHN_TO_PEDMAS_4_MISSING_NUMERIC_FACTOR(R"(\(\*)");
const	std::regex	REGEX_SHN_TO_PEDMAS_5_IMPLICIT_ADD_OPERATORS(R"((\))(\())");
const	std::regex	REGEX_SHN_TO_PEDMAS_6_PRIORITIZE_COMBINING_TERMS(R"((\([^)]*\)(?:\+\([^)]*\))+))");
const	std::regex	REGEX_SPATIAL_QTY_SIMPLIFY(R"(\(.*\))");

// set the mapped relations for GrunObject preoprties to SpatialExponentValues in here. 
// if you add additional properties to GrunObject, you need to add entries for them in here.
const std::unordered_map<std::string, SpatialExponentValue> GrunObject::propertyDimensions = {
    // Linear (1D) properties (L^1)
    {"L", SpatialExponentValue::Linear},			// Length
    {"W", SpatialExponentValue::Linear},			// Width
    {"D", SpatialExponentValue::Linear},			// Depth

    // Area (2D) properties (L^2)
    {"A", SpatialExponentValue::Area},				// Area  
    
    // Volume (3D) properties (L^3)
    {"V", SpatialExponentValue::Volume}				// Volume
    
    // All other properties default to Unitless (L^0)
};

/**
 * @brief GrunShape Constructor
 * * The m_x, m_y and m_z values are interpretted differently depending on the ShapeType.
 * * | Shape Type		| m_x			| m_y			| m_z			|
 * 	 | :---				| :---			| :---			| :---			|
 * 	 | **Rectangle**	| Length		| Width			| Depth			|
 * 	 | **Triangle**		| Base			| Perp. Height	| Depth			|
 *   | **Circle**		| Radius		| Unused (0.0)	| Depth			|
 * @param typeName sets the ShapeType. Value must be valid. (See enum class ShapeType)
 * @param m_x the GrunObject's 'x' value (Horizontal plane, direction 1)
 * @param m_y the GrunObject's 'y' value (Horizontal plane, direction 2 - adjacent to m_x/direction 1 on the same plane)
 * @param m_z the GrunObject's 'z' value (Vertical	 plane, direction 3 - adjacent to both m_x and m_y)
 * @param primaryMaterial the GrunObject's primary material - this will eventually be a string value fed in from a database somewhere
 * @param areaType sets the AreaType. Value must be valid, Horizontal or Vertical. (See enum class AreaType)
 * @throws invalid_argument - if provided ShapeType is unknown.
 * @throws invalid_argument - if x,y or z values are invalid for calculations (<= 0)
 * @throws runtime_error - if calculated aspect ratio is NaN or is infinite
 */
GrunObject::GrunObject(const std::string &shapeType,
					   const std::string &name,
                       double x,
                       double y,
					   double z,
					   const std::string &areaType,
					   const std::string &stage)
        : m_type(shapeTypeFromString(shapeType)),
		  m_name(name),
          m_x(x),
          m_y(y),
		  m_z(z),
		  m_areaType(areaTypeFromString(areaType)),
		  m_stage(stage)
{
	// class constructor
	// 1. check for member values, if some values are missing use defaults or output error?
	// 2. based on the ShapeType, can we calculate any relevant information about the object using the given member values?

	if (m_type == ShapeType::Unknown)
		throw std::invalid_argument("invalid shape m_type: " + shapeType);
	
	// calculate all relevant information possible for each ShapeType
	switch (m_type)
	{
		case ShapeType::Rectangle:
		{
			// calculate Rectangle's data
			double larger_side;			// for aspect ratio
			double smaller_side;		// for aspect ratio
			double aspectRatio;

			// calculate area of GrunObject based on m_areaType - switch for Vertical type, otherwise fallback to default Horizontal area calculation
			m_area = (m_areaType == AreaType::Vertical) ? (m_x * m_z) : (m_x * m_y);

			// calculate aspect ratio of the GrunObject based on the m_areaType
			if (m_areaType == AreaType::Vertical)
			{
				// zero check - values must be > 0
				if (m_x <= 0 || m_z <= 0) 
				{
					throw std::invalid_argument(std::format("Dimensions for aspect ratio must be positive. Values were: m_x: {}, m_z: {}", m_x, m_z));
				}
				larger_side		= std::max(m_x, m_z);
				smaller_side	= std::min(m_x, m_z);

				// paranoid, skiiered as double check for division by 0
				if (smaller_side == 0)
				{
					throw std::invalid_argument(std::format("Smaller dimension for aspect ratio calculation equals 0. Calculation not possible. Value was: smaller_side: {}", smaller_side));
				}
			}
			else if (m_areaType == AreaType::Horizontal)
			{
				// zero check - values must be > 0
				if (m_x <= 0 || m_y <= 0) 
				{
					throw std::invalid_argument(std::format("Dimensions for aspect ratio must be positive. Values were: m_x: {}, m_y: {}", m_x, m_y));
				}
				larger_side		= std::max(m_x, m_y);
				smaller_side	= std::min(m_x, m_y);

				// skiiered sanity check for division by 0
				if (smaller_side == 0)
				{
					throw std::invalid_argument(std::format("Smaller dimension for aspect ratio calculation equals 0. Calculation not possible. Value was: smaller_side: {}", smaller_side));
				}
			}
			aspectRatio = larger_side / smaller_side;

			// post calculation sanity check
			// check for NaN (Not a Number) or infinity, which can happen in edge cases.
			if (std::isnan(aspectRatio) || std::isinf(aspectRatio)) {
				throw std::runtime_error(std::format("Resulting aspect ratio is invalid (NaN or infinity). Value was: {}", aspectRatio));
			}
			
			m_aspectRatio = aspectRatio;

			// calculate volume of GrunObject - if any of these values are 0 the volume will be zero but that's for the user to deal with
			m_volume = m_x * m_y * m_z;
			break;
		}
		case ShapeType::Triangle:
		{
			// calculate Triangle's data
			// a Triangle's area is best calculated by using the Base and Perpendicular Height of the triangle. This eliminates the need to know angles and the length of all sides.
			// therefore, in a GrunObject that is ShapeType::Triangle, m_x represents the Triangle's Base, and m_y represents the Perpendicular Height (PH)
			// m_z is the object's Depth to make it, so it is a Triangular 'block', not a pyramid
			// we will also make the ShapeType::Triangle respect the AreaType

			// calculate the area, determined by the AreaType
			m_area		= (m_areaType == AreaType::Vertical) ? (m_x * m_z) / 2 	: (m_x * m_y) / 2;
			m_volume	= (m_areaType == AreaType::Vertical) ? (m_area * m_y)	: (m_area * m_z);
			break;
		}
		case ShapeType::Circle:
		{
			// calculate Circle's data 
			// The ShapeType::Circle uses m_x value as its _radius_. m_y is not used, and m_z is the circle's Depth which creates a cylinder object in 3D.
			// The Circle does not respect the AreaType value. It's area is only ever calculated on the circular face of the Object.

			// calculate the Circle's data
			// area
			m_area			= PI * m_x * m_x;
			m_circumference	= PI * (2 * m_x);
			m_volume		= m_area * m_z;
			break;
		}
		default:
			break;
	}
}

/**
 * @brief Normalizes the property name by converting it to uppercase.
 * This ensures case-insensitive lookup (e.g., "l" becomes "L", "surfacearea" becomes "SURFACEAREA").
 */
std::string normalizePropertyName(const std::string& name) {
    std::string upper = name;
    // std::transform applies the lambda (toupper) to every character in the range
    std::transform(upper.begin(), upper.end(), upper.begin(), 
        [](unsigned char c){ return std::toupper(c); });
    return upper;
}

/**
 * @brief Fetch the Spatial Exponent Value of a specific GrunObject property
 */
SpatialExponentValue GrunObject::getSpatialUnit(const std::string& propertyName) {
    // normalize the input token
    std::string normalizedName = normalizePropertyName(propertyName);
    
    auto it = propertyDimensions.find(normalizedName);
    
    if (it != propertyDimensions.end()) {
        return it->second;
    } else {
        // Warning is printed, and Unitless (0) is returned as default safety
        std::println("Warning: Unknown property '{}'. Assuming Unitless (L^0) for dimensional inference.", propertyName);
        return SpatialExponentValue::None; 
    }
}

/**
 * @brief Returns the integer value of a SpatialExponentValue value
 */
int GrunObject::asInt(SpatialExponentValue unit) {
    // Casting the strongly-typed enum value to its underlying integer type
    return static_cast<int>(unit);
}

/**
 * @brief Adds a GrunItem to the GrunObject's m_items
 * @param <name> item's name/identifier
 * @param <relationship> SHN relationship this GrunItem instance has with it's owning GrunObject
 * @param quantityFormula formula applied to result of relationship to calculate a quantity of this Item (empty by default)
 * @param units item's unit of measure ("unit(s)" by default)
 * @param primaryLabourFormula formula applied to quantity of this Item to calculate Primary Labour quantity (empty by default)
 * @return true if successful, false if failure.
 */
bool GrunObject::addGrunItem(std::string name, int libraryID, std::string relationship, std::string relComment, std::string quantityFormula, std::string units, std::string primaryLabourFormula)
{
	// zero check
	GrunItem newItem(name, libraryID, relationship, relComment, quantityFormula, units, primaryLabourFormula);
	calculateGrunItemData(newItem);
	m_items.emplace_back(newItem);
	return true;
}

/**
 * @brief addGrunItem overload
 * 
 * @param	item		(GrunItem)	- the GrunItem object to add to this GrunObject 
 * @return				(GrunItem&)	- returns a reference to the GrunItem _in its memory location the GrunObject's m_items vector_
 */
GrunItem& GrunObject::addGrunItem(GrunItem item)
{
	// check if the object we were passed has valid data or not, return out false if it's INVALID
	if (item == GrunItem::INVALID) 
	{
		std::println("Invalid GrunItem passed to GrunObject::addGrunItem(GrunItem &item)");
		return const_cast<GrunItem&>(GrunItem::INVALID);
	}

	// check the item has a library id value
	if (item._libraryId <= 0)
	{
		std::println("Library ID of GrunItem passed to GrunObject::addGrunItem(GrunItem &item) was <= 0 so returning false");
		return const_cast<GrunItem&>(GrunItem::INVALID);
	}

	this->calculateGrunItemData(item);
	m_items.push_back(item);

	// return a reference to the item we just added _in its memory location in the m_items vector_
	return m_items.back();
}

/**
 * @brief Adds a new RelationshipValues set to a referenced GrunItem
 */
size_t GrunObject::addGrunItemRelationship(GrunItem& item, const std::string &relationship, const std::string &relComment)
{
	// adds a new relationship set to the GrunItem, recalculates the GrunItem's values and returns the new total number of relationships the GrunItem has
	size_t numberofRelationships = item.addRelationshipValues(relationship, relComment);
	calculateGrunItemData(item);
	return numberofRelationships;
}

/**
 * @brief Removes an existing RelationshipValues set from the referenced GrunItem
 * 
 */
size_t GrunObject::rmGrunItemRelationship(GrunItem &item, const std::string &relationship, const std::string &relComment)
{

	return size_t();
}

/**
 * @brief	Finds index location(s) of GrunItem(s) in m_items vector that have _itemNames matching findItemName argument
 * @param	findItemName	(std::string)			- the item name to search for
 * @param	useExactSearch	(bool)					- toggles if we should match the findItemName string exactly or not (default: true). If this is set to false, the findItemName string can be used to find itemNames that CONTAIN findItemName, or * and ? wildcard chars can be used for broader searching options.
 * @return	vector of ints	(std::vector<size_t>)	- returns a vector of unsigned integers that represent any locations in the m_items vector that have matching item names
 */
std::vector<size_t> GrunObject::findGrunItemByItemName(const std::string& findItemName, bool useExactSearch) const
{
	std::vector<size_t>	indices;
	if (findItemName.empty()) { return indices; }

	std::string			pattern		= findItemName;
	std::string			regexStr	= cabji::wildcardsToRegexReady(pattern, useExactSearch);

	// find matching indices using straight string comparison (exact search)
	if (useExactSearch && pattern.find('*') == std::string::npos && pattern.find('?') == std::string::npos)
	{
		// if using exact search AND no wildcards are found, just do simple string comparison to check the matching
		for (size_t i = 0; i < m_items.size(); ++i)
		{
			if (m_items[i]._itemName == findItemName)
			{
				indices.push_back(i);
			}
		}
		return indices;
	}

	// find matching indices using broader regex search
	std::regex query(regexStr, std::regex_constants::icase);	// does case insensitive search
	for (size_t i = 0; i < m_items.size(); ++i)
	{
		if (std::regex_search(m_items[i]._itemName, query))
		{
			indices.push_back(i);
		}
	}
	return indices;
}

/**
 * @brief Search for entries in a GrunItem's _itemCoreValues vector that have a relationship matching the parameters supplied
 * @param itemIndex			(size_t)		- the index location of the GrunItem in the GrunObject's m_items vector
 * @param relationship		(std::string)	- the relationship string to search for
 * @param relComment		(std::string)	- the relationship comment to search for
 * @param useExactSearch	(bool)			- toggle using broad or exact search (default: true - uses exact searching)
 * @note finding a relationship can be based on either the relationship string, the comment string, or both as the relationship string can be empty
 */
std::vector<RelationshipSearchResult> GrunObject::findRelationshipByStrings(const size_t itemIndex, const std::string &relationship, const std::string &relComment, const bool useExactSearch) const
{
	// out of bounds check
	if (itemIndex >= m_items.size()) { return {}; }

	std::vector<RelationshipSearchResult> results;
	const std::vector<RelationshipValues>& coreValues = m_items[itemIndex]._itemCoreValues;
	
	// check if we're just looking for empty values first
	if (relationship.empty() && relComment.empty())
	{
		// there's nothing to search for, just look for any entries with empty relationship and relComment and return their indices
		for (size_t i = 0; i < coreValues.size(); ++i)
		{
			if (coreValues[i].relationship.empty() && coreValues[i].relComment.empty())
			{
				results.push_back({itemIndex, i});
			}
		}
		return results;	// return early to skip the regex processing
	}

	// find matching indices using straight string comparison (exact search)
	bool hasWildcards	= relationship.find_first_of("*?") != std::string::npos
						  || relComment.find_first_of("*?") != std::string::npos;

	// if using exact search AND no wildcards are found, just do simple string comparison to check the matching
	if (useExactSearch && !hasWildcards)
	{
		for (size_t i = 0; i < coreValues.size(); ++i)
		{
			// only match is the search term is not empty, and matches exactly
			bool relMatch = !relationship.empty()	&& coreValues[i].relationship == relationship;
			bool comMatch = !relComment.empty()		&& coreValues[i].relComment == relComment;
			if (relMatch || comMatch) {	results.push_back({itemIndex, i});}
		}
		return results;
	}

	// find matching indices using broader regex/wildcard search
	std::regex relRegex(cabji::wildcardsToRegexReady(relationship, useExactSearch), std::regex_constants::icase);				// does case insensitive search
	std::regex comRegex(cabji::wildcardsToRegexReady(relComment, useExactSearch), std::regex_constants::icase);
	for (size_t i = 0; i < coreValues.size(); ++i)
	{
		// logic: only run the regex_search if the search parameter was provided
		bool	relMatch	= !relationship.empty()	&& std::regex_search(coreValues[i].relationship, relRegex);
		bool	comMatch	= !relComment.empty()	&& std::regex_search(coreValues[i].relComment, comRegex);
		if (relMatch || comMatch)
		{
			results.push_back({itemIndex, i});
		}
	}
	return results;
}

std::vector<RelationshipSearchResult> GrunObject::findRelationshipByStrings(const std::vector<size_t> itemIndices, const std::string &relationship, const std::string &relComment, const bool useExactSearch) const
{
	// zero-check
	if (itemIndices.empty()) { return {}; }
	std::vector<RelationshipSearchResult>	allResults;
	// loop each index in itemIndices and pass the arguments to the overloaded version that does the grunt work
	for (size_t index : itemIndices)
	{
		auto	itemResults = findRelationshipByStrings(index, relationship, relComment, useExactSearch);
		std::ranges::copy(itemResults, std::back_inserter(allResults));
	}
	return allResults;
}

/**
 * @brief Calculates the aspect ratio of the GrunObject (if possible). Bases the aspect ratio calculation on the GrunObject's AreaType value.
 * @return The ratio of the longer side to the shorter side as a double (e.g., 1.5 for a 3x2 shape).
 * @throws std::invalid_argument if either dimension is zero, negative, or invalid.
 */
double GrunObject::getAspectRatio()
{ 
	return m_aspectRatio; 
}

GrunItem& GrunObject::getGrunItemByIndex(const int index)
{
	if (index >= m_items.size())
	{ 
		throw std::out_of_range("GrunObject::getGrunItemByIndex(): Index out of bounds - index = " + index);
	}
	else { return m_items[index]; }
}

/**
 * @brief Gets the GrunObject's name from member m_name
 * @return std::string containing the GrunObject's name
 */
std::string GrunObject::getObjectName()
{
	return this->m_name;
}

/**
 * @brief Get an owned GrunItem's ItemQtyTotal by index
 * 
 * @param index 		(size_t)	- the index of the GrunItem in its owning GrunObject m_item[] vector
 * @param getRounded	(bool)		- toggle if you want the unrounded or rounded total quantity
 * @return (ItemAndTotal)			- an ItemAndTotal struct instance holding the desired data
 */
ItemAndTotal GrunObject::getItemQtyTotal(const size_t& index, const bool& getRounded)
{
	// zero-check and out-of-bounds check
	if (index < 0 || index >= m_items.size()) { return ItemAndTotal(); }

	ItemAndTotal	returnObj;
	// get the GrunItem's data from m_item[index]
	returnObj.itemName	= m_items[index]._itemName;
	returnObj.itemTotal._total	= (getRounded) ? m_items[index]._itemTotalQuantityRounded : m_items[index]._itemTotalQuantity;
	return returnObj;
}

std::vector<ItemAndTotal> GrunObject::getItemQtyTotals(const bool& getRounded)
{
	// zero-check
	if (m_items.size() == 0) { return std::vector<ItemAndTotal>(); }

	std::vector<ItemAndTotal>	returnVec;
	// loop the size of the m_items vector, get the ItemQtyTotal data we want using getItemQtyTotal, push it into returnVec
	for (size_t i = 0; i < m_items.size(); i++)
	{
		returnVec.push_back(getItemQtyTotal(i, getRounded));
	}
	return returnVec;
}

double GrunObject::getObjectProperty(const std::string propertyName)
{
	// quick and dirty hack code for testing only, needs to support AreaType recognition
	if (propertyName == "length")	{ return m_x; }
	if (propertyName == "width")	{ return m_y; }
	if (propertyName == "depth")	{ return m_z; }
	if (propertyName == "area") 	{ return m_area; }
	if (propertyName == "volume")	{ return m_volume; }
	return 0.0;
}

size_t GrunObject::getTotalOfGrunItems()
{
	return m_items.size();
}

/**
 * @brief Gets the GrunObject's total cost (in cents - no decimal) to be constructed
 * 
 * @param getRounded 
 * @return long long 
 */
long long GrunObject::getTotalCostOfObject(const bool& getRounded)
{
	long long result = 0;
	for (const auto& item : m_items)
	{		
		(getRounded) ? result += item._itemTotalQuantityRounded * item._itemCostPerUnitCents : result += item._itemTotalQuantity * item._itemCostPerUnitCents;
		std::println("{} cost: {} qty: {}",item._itemName,item._itemCostPerUnitCents,item._itemTotalQuantity);
	}
	std::println("result cost: {}", result);
	return result;
}

/**
 * @brief Removes one (first found) or all GrunItems with the specified name from the GrunObject's m_items.
 * @param itemName The name of the GrunItem(s) to remove.
 * @param removeAll If true, removes all items with the name; otherwise, removes only the first match found using iterative search.
 * @return The number of items removed.
 */
size_t GrunObject::removeGrunItem(const std::string& itemName, bool removeAll)
{
    size_t initialSize = m_items.size();

    // create a lambda predicate to check if the itemName matches the _itemName member of the worked on GrunItem
    auto itemNameIsFound = [itemName](const GrunItem& item) {
        return item._itemName == itemName;
    };

	// if we want to remove all instances of GrunItems with _itemName matching itemName
    if (removeAll)
    {
        // std::remove_if moves all elements that match the predicate to the end, 
        // and returns an iterator to the location of the first element to 'remove'
        auto it = std::remove_if(m_items.begin(), m_items.end(), itemNameIsFound);
        
        // std::vector::erase actually removes the elements located from the iterator 'it' to the end of the vector.
        m_items.erase(it, m_items.end());

        // Return the count of removed items
        return initialSize - m_items.size();
    }
    else // Remove only the first match
    {
        // 1. Find the first item that matches the predicate
        auto it = std::find_if(m_items.begin(), m_items.end(), itemNameIsFound);

        if (it != m_items.end())
        {
            // 2. Erase that single element
            m_items.erase(it);
            return 1;
        }
        return 0; // No item found
    }
}

/**
 * @brief Removes a single GrunItem at a specific index.
 * @param index The zero-based (integer value) index of the item to remove.
 * @return true if the item was successfully removed, false otherwise (e.g., index out of bounds).
 */
bool GrunObject::removeGrunItem(size_t index)
{
    if (index < m_items.size())
    {
        // std::vector::erase is efficient for a single element removal
        m_items.erase(m_items.begin() + index);
        return true;
    }
    return false; // Index was out of bounds
}

/**
 * @brief Removes multiple GrunItems at specific indices.
 * @param indices A std::vector of zero-based (integer value) indices of the items to remove.
 * @return the number of elements successfully removed by the function
 */
size_t GrunObject::removeGrunItem(std::vector<size_t> indices)
{
	// dev-note: this function removes multiple GrunItems in m_items as long as the item's index value is valid.
	// it removes the items from m_items in reverse order.

	// zero check
	if (indices.size() == 0) { return 0; }
	size_t returnVal = 0;
	// sort vector values in descending order
	std::sort(indices.begin(),indices.end(), std::greater<int>());
	for (auto index : indices)
	{
		if (index < m_items.size())
		{
			m_items.erase(m_items.begin() + index);
			returnVal++;
		}
	}
	return returnVal;
}

/**
 * @brief Returns the number of items in the m_items vector
 * 
 * @return size_t 
 */
size_t GrunObject::size()
{
	return m_items.size();
}

/**
 * @brief Gets a string that will describes the items in the GrunObject (mostly for debugging output)
 * @param	dateFormat		(std::string)	- a str that defines the format for the date used for the LKGW calculated time for GrunItems
 * @param	itemNameWidth	(size_t)		- an unsigned
 * @return	std::string showing information about GrunItems in the GrunObject
 */
std::string GrunObject::getGrunItemListInfoAsString(const std::string dateFormat, const size_t itemNameWidth)
{
	std::string returnVal;
	if (!m_items.empty()) { returnVal = "GrunItems in Element '" + m_name +"'\n"; }
	for (auto item : m_items)
	{
		returnVal += std::format("{:<{}} ",item._itemName.substr(0,itemNameWidth), itemNameWidth);
		
		// deal with item's CoreValue vector
		int			numRels			= item.getNumberOfRelationships();
		std::string	osRelationship	= "";
		std::string osSpatialQty	= "";
		std::string osSpatialUnit	= "";
		std::string osItemQty		= "";

		for (const auto coreValueSet : item._itemCoreValues)
		{
			osRelationship	+= coreValueSet.relationship											+ ", ";
			osSpatialUnit	+= spatialExponentValueToString(coreValueSet.spatialUnit).substr(0,9)	+ ", ";
		}

		osSpatialQty	= std::format("{:7.2f} ",item._spatialTotalQuantity);
		osItemQty		= std::format("{:7.2f} ",item._itemTotalQuantity);
		// remove trailing ", "
		osRelationship	= osRelationship.substr(0, (osRelationship.length() - 2));
		osSpatialUnit	= osSpatialUnit.substr(0, (osSpatialUnit.length() - 2));

		returnVal += std::format("Rel: {:<12} ",osRelationship.substr(0,12)) + 			 
					 std::format("SQ: {}",osSpatialQty) + 
					 std::format("SU:{:<11} ",osSpatialUnit.substr(0,11)) + 
					 std::format("{:<10}","Item.Qty: ") + 
					 std::format("{:<9} ",osItemQty) +
					 std::format("{:<8} ",item._itemQuantityUnits.substr(0,8)) +
					 std::format("(R:{:>7.2f}) ",item._itemTotalQuantityRounded) +
					 std::format("{:<10} ","P.Labour:") + 
					 std::format("{:>6.2f} ", item._itemTotalPrimaryLabour) + 
					 std::format("{:<6} "," hr(s)") +
					 std::format("{:<6} ","L.CT:") + 
					 std::format("{:>}", item.getCalculatedTimeString(item._itemLKGWCalculated, dateFormat).substr(0,13)) +
					 "\n";
	}
	return returnVal;
}

/**
 * @brief Formats the calculated GrunObjectTotals into a readable string.
 * @return A string containing all aggregated totals.
 */
std::string GrunObject::getGrunObjectTotalsInfoAsString() const
{
	std::stringstream ss;
	
	// Helper array to hold the display names for the map categories
	const std::string categories[] = {
		"Primary Labour Total",
		"Spatial Material Totals (By Spatial Unit)",
		"Item Unit Totals (By Item Name)"
	};
	
	// Access all three maps via the GrunObjectTotals instance
	const auto& totals = m_objectTotals;
	
	// Array of pointers to the three map members for iteration
	auto totalMaps = GrunObjectTotals::TOTALS_PTRS;

	ss << std::format("\n--- {} Aggregated Totals ---\n", m_name);
	
	for (size_t i = 0; i < GrunObjectTotals::getMapCount(); ++i)
	{
		// 1. Get the current map (read-only reference)
		const auto& currentMap = totals.*totalMaps[i];

		ss << std::format("\n  [ {} ]\n", categories[i]);
		
		if (currentMap.empty())
		{
			ss << std::format("    (No aggregated data)\n");
			continue;
		}

		// 2. Iterate through the key-value pairs in the map
		for (const auto& [key, totalAndUnit] : currentMap)
		{
			// Output format: [Key] = Total (Unit)
			// Using {:.<35} for left-justified key with padding of dots
			ss << std::format("    {:.<35} = {}\n", key, totalAndUnit);
		}
	}
	
	return ss.str();
}

/**
 * @brief combines a value (lhs) with an additional math formula (formula)
 */
double GrunObject::applyFormula(double lhs, const std::string &formula)
{
	if (formula.empty()) return lhs;

	std::regex formula_pattern("\\s*([\\+\\-\\*\\/])\\s*(\\d*\\.?\\d+)");
	std::smatch match;

	if (std::regex_match(formula, match, formula_pattern))
	{
		char op = match[1].str()[0];
		double val = std::stod(match[2].str());
		switch (op)
		{
			case '+': return lhs + val;
			case '-': return lhs - val;
			case '*': return lhs * val;
			case '/': return (val != 0) ? lhs / val : lhs;
		}
	}
	return lhs;
}

/**
 * @brief Calculates a GrunItem's member values
 * @param item (GrunItem) - the GrunItem to calculate the data for
 * @return true if successful, false if failure occurs
 */
bool GrunObject::calculateGrunItemData(GrunItem &item)
{
	// first interpret the item's relationship to calculate Spatial Values and Item Quantity 
	interpretGrunItemSpatialValues(item);
	interpretGrunItemItemQuantity(item);

	// loop the item._itemCoreValues vector and total up the calculated values
	/*
	_spatialTotalQuantity
	_itemTotalQuantity
	_itemTotalPrimaryLabour
	*/

	// reset totals to 0 otherwise we end up with incorrect totals for items with multiple relationships
	item._spatialTotalQuantity	= 0;
	item._itemTotalQuantity		= 0;
	for (const auto coreValueSet : item._itemCoreValues)
	{
		item._spatialTotalQuantity	+= coreValueSet.spatialQuantity;
		item._itemTotalQuantity		+= coreValueSet.itemQuantity;
	}

	// next, calculate all the GrunItem's member values
	item._itemTotalPrimaryLabour		= applyFormula(item._itemTotalQuantity,item._itemPrimaryLabourFormula);
	item._itemTotalQuantityRounded		= cabji::roundToStep(item._itemTotalQuantity,item._itemRoundUpFactor);
	item._itemWasteAllowance			= item._itemTotalQuantity * item._itemWasteFactor;
	item._itemItemizedProfit			= item._itemTotalQuantity * item._itemItemizedProfitFactor;
	item._itemLKGWCalculated			= std::chrono::system_clock::now();
	return true;
}

int GrunObject::calculateGrunObjectTotals() {
    for (int i = 0; i < GrunObjectTotals::getMapCount(); i++) 
	{
        auto memberPtr = GrunObjectTotals::TOTALS_PTRS[i];
        auto& currentMap = m_objectTotals.*memberPtr;

        for (const auto& item : m_items) 
		{
            // Case 0 and Case 2 only care about the ITEM as a whole
			// this if conditional is where we set what item members to look at for each totalling aggregation
            if (i == 0 || i == 2) 
			{
                std::string aggregationKey = (i == 0) ? m_name : item._itemName;
                double itemQuantity = (i == 0) ? item._itemTotalPrimaryLabour : item._itemTotalQuantityRounded;
                std::string itemUnit = (i == 0) ? item._itemPrimaryLabourUnits : item._itemQuantityUnits;

                if (itemQuantity != 0.0 && !aggregationKey.empty()) 
				{
                    TotalAndUnit& entry = currentMap[aggregationKey];
                    entry._total += itemQuantity;
                    if (entry._unit == "unit(s)" || entry._unit.empty()) entry._unit = itemUnit;
                }
            }
            // Case 1: Spatial Totals (Iterate through the new relationship vector)
            else if (i == 1) {
                for (const auto& rel : item._itemCoreValues) 
				{
					// Case 1 distinguish between different spatial unit types "Timber (m)" and "Timber (m3)"
					std::string aggregationKey = item._itemName + " (" + spatialUnitToSuffix(rel.spatialUnit) + ")";
                    double relQuantity = rel.spatialQuantity;

                    if (relQuantity != 0.0) {
                        TotalAndUnit& entry = currentMap[aggregationKey];
                        entry._total += relQuantity;
                        
                        // Set the unit based on the spatial dimension (m, m2, m3, etc)
                        if (entry._unit == "unit(s)" || entry._unit.empty()) {
                            entry._unit = spatialUnitToSuffix(rel.spatialUnit); 
                        }
                    }
                }
            }
        }
    }
    return 0;
}

/**
 * @brief Calculates the totals for this Grunobject
 * 
 * @return true 
 * @return false 
 */
bool GrunObject::calculateTotals()
{
	m_TOTAL_COST_ItemsCalculated		= 0;
	m_TOTAL_COST_ItemsRounded			= 0;
	m_TOTAL_HOURS_LabourCalculated		= 0;
	
	// zero-check
	if (m_items.size() > 0) 
	{
		// loop the items vector
		for (const auto& item : m_items)
		{
			m_TOTAL_COST_ItemsCalculated	+= item._itemTotalQuantity 			* item._itemCostPerUnitCents;
			m_TOTAL_COST_ItemsRounded 		+= item._itemTotalQuantityRounded	* item._itemCostPerUnitCents;
			m_TOTAL_HOURS_LabourCalculated 	+= item._itemTotalPrimaryLabour;
			// dev-note: labour TOTAL_COST is not calculable at the GrunObject scope because the hourly rate is not available here
			// dev-note: rounded labour total is not calculable at the GrunObject scope because it will eventually require 'Sessioning' and sessioning is not implemented in the MVP
		}
	}
	
	return true;
}

// looks at a GrunItem's unit of measure and tries to determine which SpatialExponentValue it should be
SpatialExponentValue GrunObject::mapUnitToSpatialExponent(const std::string& unit) const
{
    std::string lowerUnit = unit;
    std::transform(lowerUnit.begin(), lowerUnit.end(), lowerUnit.begin(),
        [](unsigned char c){ return std::tolower(c); });

    if (lowerUnit.find("m3") != std::string::npos || lowerUnit.find("cubic") != std::string::npos)
    {
        return SpatialExponentValue::Volume;
    }
    else if (lowerUnit.find("m2") != std::string::npos || lowerUnit.find("square") != std::string::npos)
    {
        return SpatialExponentValue::Area;
    }
    else if (lowerUnit.find("m") != std::string::npos || lowerUnit.find("lineal") != std::string::npos)
    {
        return SpatialExponentValue::Linear;
    }
    
    return SpatialExponentValue::None;
}

/**
 * @brief This function will replace Grun Object Token characters (LWDAV etc) with the appropriate numeric value from the owning Grun Object's property
 * @param (std::string) - the SHN relationship string to parse
 * @return (std::string) - returns the parsed string with tokens replaced with numbers
 */
std::string GrunObject::substituteRelationshipTokens(const std::string& relationship) const
{
    std::string substituted_relationship = relationship;

    // 1. Define the substitution values (ensure they are in string format for replacement)
    // We prioritize using the explicit object dimension members for L, W, D.
    std::string L_val = std::to_string(m_x);
	std::string W_val = std::to_string(m_y); // Use m_y for W
	std::string D_val = std::to_string(m_z); // Use m_z for D
	if (m_areaType == AreaType::Vertical)
	{
		// the GrunObject is a "vertical" area (like a wall)
		W_val = std::to_string(m_z); // Use m_z for W
		D_val = std::to_string(m_y); // Use m_y for D
	}
	
    std::string A_val = std::to_string(m_area);
    std::string V_val = std::to_string(m_volume);
    std::string P_val = std::to_string(2.0 * (m_x + m_y)); // Assuming a rectangular perimeter for 'P'
    std::string C_val = std::to_string(m_circumference);

    // 2. Perform substitutions using regex (Need to be careful about order, V and A first)
    // The substitution logic must ensure that "2L1W" becomes "2*L + 1*W" where no operator is present
    
    // Pattern to find "Number[TOKEN]" or just "[TOKEN]" (e.g., 2L, 1W, A)
    // This is complex, but the safest way is to substitute tokens first, then insert operators.
    
    // Simple Substitution: Replace V, A, L, W, D, P, C tokens with their values.
    // Replace complex tokens first (V, A), then linear (L, W, D).

    std::vector<std::pair<std::string, std::string>> replacements = {
        {"V", V_val},
        {"A", A_val},
        {"P", P_val},
        {"C", C_val},
        // The L, W, D replacements need to handle "2L" vs "L" vs "L*W"
        {"D", D_val},
        {"W", W_val},
        {"L", L_val}
    };
    
    // Perform substitution for tokens followed by non-number/non-operator characters (e.g., L)
    for (const auto& rep : replacements)
    {
        // Use a regex word boundary \b to ensure 'L' isn't replaced inside 'ROLLS' (if that were a token)
        // Since your tokens are single capital letters, a simple search-and-replace is often adequate
        size_t pos = substituted_relationship.find(rep.first);
        while (pos != std::string::npos)
        {
            substituted_relationship.replace(pos, rep.first.length(), rep.second);
            pos = substituted_relationship.find(rep.first, pos + rep.second.length());
        }
    }

    return substituted_relationship;
}

/**
 * @brief This function is responsible for all operations involved in determining all values associated with a GrunItem's Spatial Values. It will look at the item's _relationship value and determine all necessary Spatial related values and GrunItem member values when this function is called.
 * @param (GrunItem) - the Grun Item instance to interpret/calculate the Spatial Values for
 * @return (bool) - false if unsuccessful (this will occur if the Grun Item's relationship string is faulty), true if successful.
 */
bool GrunObject::interpretGrunItemSpatialValues(GrunItem &item)
{
	bool foundARelationship = false;

// 0: loop the _itemCoreValues vector for access to relationship string(s)
	// loop the _itemCoreValues vector and handle each relationship the item has
	for (RelationshipValues& coreValueSet : item._itemCoreValues)
	{
		// data acquisition
		std::string	relationship	= coreValueSet.relationship;
		std::string	baseExpr		= "";
		// get the base expression (anything before the last occurence of an @ char, or the whole string)
		auto		atPos			= relationship.find_last_of('@');
		baseExpr = relationship.substr(0, atPos);
		std::string	resultPattern	= "";
		// use regex_replace to find significant characters in the baseExpr and put them in a string														// the replacement pattern used for a regex_replace call
		std::string	saneBaseExpr	= std::regex_replace(baseExpr, REGEX_GI_BASEEXPR_SIG_TOKENS_AND_OPS, resultPattern);
		coreValueSet.baseExpression = baseExpr;
		auto		current			= std::source_location::current();							// for debugging output if needed
		int			spatialAnchor	= 0;

		// zero-check: if relationship is empty, 'continue' to skip for loop
		if (saneBaseExpr.empty())
			continue;	// end of zero-check

		foundARelationship = true;
		
	// 1: sanitize the base expression for Spatial Unit and Spatial Value calculation

		// we are assuming the saneBaseExpr has *something* in it at this point - probably dangerous, but we did do an .empty() check
		// strip * or / from front and back of string
		while ((saneBaseExpr.front() == '*' || saneBaseExpr.front() == '/'))
		{
			// chomp the first char if its * or /
			saneBaseExpr.erase(0,1);
		}
		while ((saneBaseExpr.back() == '*' || saneBaseExpr.back() == '/'))
		{
			// chomp the last char if its * or /
			saneBaseExpr.pop_back();
		}

		// assign saneBaseExpr to the appropriate member in the item
		coreValueSet.baseExpressionIntrpForSU = saneBaseExpr;

	// 2: set a few values that are mostly only useful for bug tracking if you get bugs appearing
		
		// the saneBaseExpr with GrunObject Tokens switched out for their numeric worth
		std::string numericExpr = saneBaseExpr;

		// dev-note: we are assuming saneBaseExpr is a SANITIZED string. If you get unexpected behaviour, you should probably check the value of saneBaseExpr
		// convert GrunObject Tokens in numericExpr to their Spatial Exponent Values and convert * operators to +
		for (char& c : numericExpr) 
		{
			// convert * operators to +, else convert GrunObject Token to its SpatialExponentValue and cast the number back as a char in the numericExpr string
			if (c == '*')
				c = '+';
			else
				c = static_cast<char>(asInt(getTokenExponent(c)) + 48);	// dev-note: + 48 to correctly cast the returned int value BACK to a char
		}

		coreValueSet.baseExpressionIntrpNumeric = numericExpr;

		// get the largest number in the numericExpr, clamp to min 0 and max 3
		auto digits	= numericExpr 
					| std::views::filter(::isdigit)
					| std::views::transform([](char c) { return c - '0'; });

		if (!digits.empty())
			spatialAnchor 			= std::clamp(std::ranges::max(digits),0,3);
		coreValueSet.spatialAnchor	= static_cast<SpatialExponentValue>(spatialAnchor);						  

		// process the operators - loop through the numericExpr char by char with the index value available
		std::string	numericExprResult;
		// used to skip operators during processing so we only add numbers together
		int			skipCount = 0;
		for (auto [i, c] : std::views::enumerate(numericExpr))
		{
			// skip this iteration of the loop if skip count is counting down
			if (skipCount > 0) { skipCount--; continue; }

			std::optional<int> r;
			if (c == '+' && i > 0 && (i+1) < numericExpr.length())				// check we still have 1 more char in the numericExpr string
			{
				// if we detect a ++ operator, skip the next 2 iterations
				if ((i+2) < numericExpr.length() && numericExpr[i+1] == '+') { skipCount = 2; continue; }

				// standard logic for single '+'
				const int	lhs		= numericExpr[i-1] - '0';					// dev-note: we must use - '0' here to convert the CHAR 1 into type int 1
				const int	rhs		= numericExpr[i+1] - '0';
				r	= lhs + rhs;
			}
			if (r.has_value())
				numericExprResult += std::to_string(r.value());
			else
			{
				const bool	plusOnLeft 	= (i > 0 && numericExpr[i-1] == '+');
				const bool	plusOnRight	= (i+1 < numericExpr.length() && numericExpr[i+1] == '+');
				if (!plusOnLeft && !plusOnRight && c != '+')
					numericExprResult += c;
			}
		}

		int spatialUnit = 0;
		auto digitsAfter	= numericExprResult
							| std::views::transform([](char c) { return c - '0'; });
		if (!digitsAfter.empty())
			spatialUnit = std::ranges::max(digitsAfter);
		
		coreValueSet.spatialUnit = static_cast<SpatialExponentValue>(spatialUnit);

	// 3. calculate the Spatial Quantity by evaluating the GrunItem's base expression
		coreValueSet.spatialQuantityFormula = convertSpatialQuantitySHNToPEDMAS(coreValueSet.baseExpression);
		coreValueSet.spatialQuantityFormula = substituteRelationshipTokens(coreValueSet.spatialQuantityFormula);
		coreValueSet.spatialQuantity 		= evaluateArithmetic(coreValueSet.spatialQuantityFormula);

		// debug output
		// std::println("Debug Output in: {}",current.function_name());
		// std::print("item.rel: {:>15} ",item._relationship.substr(0,25));
		// std::print("baseExpr: {:>10} ",item._baseExpression.substr(0,20));
		// std::print("bEForSV: {:>5} ",item._baseExpressionIntprForSU.substr(0,14));
		// std::print("bEIntNum: {:>5} ",item._baseExpressionIntprNumeric.substr(0,15));
		// std::print("S.A.: {:>7} ",spatialExponentValueToString(item._spatialAnchor));
		// std::print("numericExprRes: {:>6} ",numericExprResult.substr(0,22));
		// std::print("S.U.: {:>7} ",spatialExponentValueToString(item._spatialUnit));
		// std::println("SQF: {:>25} ",item._spatialQuantityFormula.substr(0,30));
	}

	// return true to indicate interpretation was a success
	if (foundARelationship) 
		return true;
	else
		return false;
}

/**
 * @brief Parses a GrunItem's SHN to (hopefully) have correctly interpretted operators and parentheses for PEDMAS arithmetic evaluation
 * @param (std::string) - the shorthand notation as a string to parse and convert
 * @return (std::string) - the converted result string
 */
std::string GrunObject::convertSpatialQuantitySHNToPEDMAS(const std::string &shn)
{
	// dev-note: when we enter this method, shn should be the "base expression" which is anything BEFORE an @ sign in the relationship string of an item
	// therefore, argument shn should look something like: 2L1W, 3.0*2L or 5L etc.

	// zero-check
	if (shn.empty()) return std::string();

	// use regexes to process the shn string in our custom order of precedence to build the formula PEDMAS how we need it
	// dev-note: regex patterns used here are defined in constants at hte top of this file.

	std::string numericForm = shn;
	std::erase_if(numericForm, [](char c) { return std::isspace(static_cast<unsigned char>(c)); });				// strip whitespace
	numericForm	= std::regex_replace(numericForm, REGEX_SHN_TO_PEDMAS_2_NUM_FACTOR_AND_GO_TOKEN, "($1*$2)");	// puts implicit multiplication operators in
	numericForm = std::regex_replace(numericForm, REGEX_SHN_TO_PEDMAS_4_MISSING_NUMERIC_FACTOR, "(1*");			// puts "1*" in where required
	numericForm = std::regex_replace(numericForm, REGEX_SHN_TO_PEDMAS_5_IMPLICIT_ADD_OPERATORS, ")+(");			// puts implicit addition operators in

	std::smatch match;
	if (std::regex_search(numericForm, match, REGEX_SPATIAL_QTY_SIMPLIFY))
	{
		numericForm = match[0].str();
	}
	else
	{
		// dev-note: if we get here this is a problem!
		numericForm = "";
	}
	
	return numericForm;
}

bool GrunObject::interpretGrunItemItemQuantity(GrunItem &item)
{
	// zero-check
	if (item._itemCoreValues.empty()) return false;

	for (RelationshipValues& coreValueSet : item._itemCoreValues)
	{
		std::string	relationship	= coreValueSet.relationship;
		std::erase_if(relationship, [](char c) { return std::isspace(static_cast<unsigned char>(c)); });			// strip whitespace
		//relationship	= std::regex_replace(relationship, REGEX_SHN_TO_PEDMAS_0_WRAP_ALL_IN_PARENTHESES, "($1)");
		relationship	= std::regex_replace(relationship, REGEX_SHN_TO_PEDMAS_1_EXPLICIT_COMBINE_OPERATOR, ")$1(");
		relationship	= std::regex_replace(relationship, REGEX_SHN_TO_PEDMAS_2_NUM_FACTOR_AND_GO_TOKEN, "($1*$2)");
		relationship	= std::regex_replace(relationship, REGEX_SHN_TO_PEDMAS_3_AT_OPERATOR, "/$2+1");
		relationship	= std::regex_replace(relationship, REGEX_SHN_TO_PEDMAS_4_MISSING_NUMERIC_FACTOR, "(1*");
		relationship	= std::regex_replace(relationship, REGEX_SHN_TO_PEDMAS_5_IMPLICIT_ADD_OPERATORS, ")+(");
		relationship	= std::regex_replace(relationship, REGEX_SHN_TO_PEDMAS_6_PRIORITIZE_COMBINING_TERMS, "($1)");

		coreValueSet.interprettedRelationship = substituteRelationshipTokens(relationship);
		double	itemQty	= evaluateArithmetic(coreValueSet.interprettedRelationship);
		// itemQty at this point is only preliminary. check if the item has a value for its _itemQuantityFormula and apply it if it does
		if (!item._itemQuantityFormula.empty())
		{
			std::string	tempForm	= std::to_string(itemQty) + item._itemQuantityFormula;
			itemQty					= evaluateArithmetic(tempForm);
		}
		
		coreValueSet.itemQuantity = itemQty;
		// std::println("relationship: {}",relationship);
	}
	return false;
}

/**
 * @brief Evaluates an expression string into a numeric double value, obeying PEDMAS order of precedence
 * @param expression - (std::string) the expression to evaluate
 * @return (double) - the evaluated result of expression
 * */
double GrunObject::evaluateArithmetic(std::string expression)
{
	if (expression.empty()) return 0.0;

	// handle parentheses recursively
	size_t openBracket = expression.find_last_of('(');
	while (openBracket != std::string::npos)
	{
		size_t closeBracket = expression.find(')', openBracket);
		
		// safety break
		if (closeBracket == std::string::npos) 
			break;
		
		// extract the inside expression, evaluate it, and swap it back into the string
		std::string inside = expression.substr(openBracket + 1, closeBracket - openBracket - 1);
		expression.replace(openBracket, closeBracket - openBracket + 1, std::to_string(evaluateArithmetic(inside)));
		openBracket = expression.find_last_of('(');
	}

	// lambda function in GrunObject::evaluateArithmetic() that evaluates a math expression string down to its numeric result obeying the PEDMAS order of precedence
    auto performPass = [&](const std::string& opPattern) {
        std::regex pattern(R"(\(*([-+]?\d*\.?\d+)\)*\s*()" + opPattern + R"()\s*\(*([-+]?\d*\.?\d+)\)*)");
        std::smatch match;
        while (std::regex_search(expression, match, pattern)) {
            double left = std::stod(match[1].str());
            char op = match[2].str()[0];
            double right = std::stod(match[3].str());
            double res = 0.0;

            if (op == '*') res = left * right;
            else if (op == '/') res = (right != 0) ? left / right : 0.0;
            else if (op == '+') res = left + right;
            else if (op == '-') res = left - right;

            expression.replace(match.position(), match.length(), std::to_string(res));
        }
    };

    performPass(R"([\*/])"); // Pass 1: Multiplication and Division
    performPass(R"([\+\-])"); // Pass 2: Addition and Subtraction

	// final cleanup - remove any stray whitespace or parentheses before stod() is called
	expression.erase(std::remove(expression.begin(), expression.end(), '('), expression.end());
	expression.erase(std::remove(expression.begin(), expression.end(), ')'), expression.end());
    try 
	{ 
		return std::stod(expression);
	}
    catch (...)	
	{ 
		return 0.0;
	}
}

/**
 * @brief takes a string (assumed the string represents a GrunObject Token) and returns the SpatialExponentValue of the Token
 * @param (std::string_view / string) - the token (as a string type)
 * @return (SpatialExponentValue) - the SpatialExponentValue that the token was found to match best to, if any (returns 'None' if token is not known)
 */
SpatialExponentValue GrunObject::getTokenExponent(std::string_view token)
{

	if (GO_LINEAL_TOKENS.contains(token))	return SpatialExponentValue::Linear;
	if (GO_AREA_TOKENS.contains(token))		return SpatialExponentValue::Area;
	if (GO_VOLUME_TOKENS.contains(token))	return SpatialExponentValue::Volume;

    // if (token == "V") return SpatialExponentValue::Volume;
    // if (token == "A") return SpatialExponentValue::Area;
    // if (token == "L"	|| token == "W"	|| token == "D" || 
        // token == "PH"	|| token == "C"	|| token == "R") return SpatialExponentValue::Linear;   
    return SpatialExponentValue::None;
}

/**
 * @brief overloaded function to accept token as a char type
 * @param (char) - the token (as a char type)
 * @return (SpatialExponentValue) - the SpatialExponentValue that the token was found to match best to, if any (returns 'None' if token is not known)* 
 */
SpatialExponentValue GrunObject::getTokenExponent(char token)
{
	std::string s(1,token);
	return getTokenExponent(s);
}

std::string GrunObject::spatialUnitToString(SpatialExponentValue unit) {
    switch(unit) {
        case SpatialExponentValue::Linear: return "Linear Totals";
        case SpatialExponentValue::Area: return "Area Totals";
        case SpatialExponentValue::Volume: return "Volume Totals";
        default: return "Unitless Totals";
    }
}

std::string GrunObject::spatialUnitToSuffix(SpatialExponentValue unit) {
    switch(unit) {
        case SpatialExponentValue::Linear: return "m";
        case SpatialExponentValue::Area: return "m2";
        case SpatialExponentValue::Volume: return "m3";
        default: return "unit(s)";
    }
}