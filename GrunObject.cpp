#include <algorithm>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <optional>
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
GrunObject::GrunObject(const std::string &typeName,
					   const std::string &name,
                       double x,
                       double y,
					   double z,
					   const std::string &areaType,
					   const std::string &stage)
        : m_type(shapeTypeFromString(typeName)),
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
		throw std::invalid_argument("invalid shape m_type: " + typeName);
	
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
bool GrunObject::addGrunItem(std::string name, std::string relationship, std::string quantityFormula, std::string units, std::string primaryLabourFormula)
{
	// zero check
	GrunItem newItem(name, relationship, quantityFormula, units, primaryLabourFormula);
	interpretGrunItemSpatialValues(newItem);
	interpretGrunItemItemQuantity(newItem);
	calculateGrunItemData(newItem);
	m_items.emplace_back(newItem);
	return true;
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

/**
 * @brief Gets the GrunObject's name from member m_name
 * @return std::string containing the GrunObject's name
 */
std::string GrunObject::getObjectName()
{
	return this->m_name;
}

double GrunObject::getObjectProperty(const std::string propertyName)
{
	// quick and dirty hack code for testing only, needs to support AreaType recognition
	if (propertyName == "length")	{ return m_x; }
	if (propertyName == "width")	{ return m_y; }
	if (propertyName == "depth")	{ return m_z; }
	if (propertyName == "area") 	{ return m_area; }
	return 0.0;
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
 * @brief Gets a string that will describe the secondary materials in the GrunObject (mostly for debugging output)
 * @return std::string showing information about GrunItems in the GrunObject
 */
std::string GrunObject::getGrunItemListInfoAsString(const std::string dateFormat)
{
	std::string returnVal;
	if (!m_items.empty()) { returnVal = "GrunItems in Element '" + m_name +"'\n"; }
	for (auto item : m_items)
	{
		returnVal += std::format("{:<16} ",item._itemName.substr(0,16)) +
					 std::format("{:<17} ","Rel: " + item._relationship.substr(0,12)) + 
					 std::format("{:<13} ","SA: " + spatialExponentValueToString(item._spatialAnchor).substr(0,9)) + 
					 std::format("{:<9}","Spa.Qty: ") + 
					 std::format("{:>7.2f} ",item._spatialQuantity) +
					 std::format("{:<10}","Item.Qty: ") + 
					 std::format("{:>7.2f} ",item._itemQuantity) +
					 std::format("{:<8} ",item._itemQuantityUnits.substr(0,8)) +
					 std::format("(R:{:>7.2f}) ",item._itemQuantityRounded) +
					 std::format("{:<10} ","P.Labour:") + 
					 std::format("{:>6.2f} ", item._itemPrimaryLabour) + 
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
		"Spatial Material Totals (By Relationship)",
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

bool GrunObject::calculateGrunItemData(GrunItem &item)
{
	// when this function runs, we should already have the GrunItem's spatial values and item quantity
	// calculate numerous GrunItem attributes from the existing values

	// zero-check: even if itemQuantity is 0, we still need to set everything else to 0.
	item._itemPrimaryLabour			= applyFormula(item._itemQuantity,item._itemPrimaryLabourFormula);
	item._itemQuantityRounded		= cabji::roundToStep(item._itemQuantity,item._itemRoundUpFactor);
	item._itemLKGWCalculated		= std::chrono::system_clock::now();
	return true;
}

int GrunObject::calculateGrunObjectTotals()
{
	// loop through the pointer-to-members in m_totalsPtrs using the number of members assigned in the GrunObjectTotals struct
	for (int i = 0; i < GrunObjectTotals::getMapCount(); i++)
	{
		// define the destination for our data
		// set which member of GrunObjectTotals we are pointing to using the value of i as the index of the GrunObjectTotals::TOTAL_PTRS[] array
		auto		memberPtr	= GrunObjectTotals::TOTALS_PTRS[i];
		auto&		currentMap	= m_objectTotals.*memberPtr;

		for (const auto& item : m_items)
		{
			// Temporary variables to hold the aggregation data determined by the switch
			std::string aggregationKey	= "";
			double		itemQuantity	= 0.0;
			std::string	itemUnit		= "";

			switch (i)
			{
				case 0:
					// Case 0 (Labour Total): Aggregate everything into a single grand total, keyed by the object's name.
					aggregationKey	= m_name;
					itemQuantity	= item._itemPrimaryLabour;
					itemUnit		= item._itemPrimaryLabourUnits;
					break;
				
				case 1:
					// Case 1 (Spatial Totals): Aggregate by the relationship unit (e.g., 'A', '2W', 'V')
					aggregationKey	= item._relationship;
					itemQuantity	= item._spatialQuantity;
					itemUnit		= item._relationship; // Use relationship string as the unit/identifier
					break;

				case 2:
					// Case 2 (Item Unit Totals): Aggregate by the Item Name.
					aggregationKey	= item._itemName; 
					itemQuantity	= item._itemQuantity;
					itemUnit		= item._itemQuantityUnits; // Use the base item unit
					break;
				
				default:
					break;
			}
			
			// 2. Aggregation step: Use the determined key to look up (or create) the entry in the map
			if (itemQuantity != 0.0 && !aggregationKey.empty()) 
			{
				// Get a reference to the specific TotalAndUnit structure we want to update.
				// operator[] will create a new TotalAndUnit entry if the key doesn't exist.
				TotalAndUnit& entry = currentMap[aggregationKey];

				// Aggregate the total value
				entry._total += itemQuantity;

				// Set the unit: Only set it if the entry is currently using the default/empty unit.
				// This ensures the unit from the first item aggregated is used for the key.
				if (entry._unit == "unit(s)" || entry._unit.empty()) {
					entry._unit = itemUnit;
				}
			}
		}
	}
	return 0;
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

// decides where implied multiplication operators should be inserted into a segement/term from a relationship string
std::string GrunObject::injectImplicitOperators(std::string &segment) 
{
	// data acquisition
	// token_start holds the first GrunObject Token in the segment (term) if one is found
    auto token_start = std::find_if(segment.begin(), segment.end(), [](char c) {
        return std::isupper(static_cast<unsigned char>(c));
    });

    // zero-check
	// if no GrunObject Token is found, (the segment is a direct quantity like "8" or an explicit operator with a trailing value like "@1.2") return it, as is
    if (token_start == segment.end()) return segment;

    // separate prefix and token
    std::string prefix(segment.begin(), token_start);
    std::string token(token_start, segment.end());

    // handle empty prefix (e.g., "L" -> "1*L")
    if (prefix.empty()) {
        return "1*" + token;
    }

    // check if we need to inject '*'
    char lastChar = prefix.back();
    
    // if lastChar is a digit or decimal, inject '*'
    if (std::isdigit(static_cast<unsigned char>(lastChar)) || lastChar == '.') {
        return prefix + "*" + token;
    }

    // if it ends in an anything else (assumed: /, +, -, *), join the pieces back together, example: "1.5/L" -> "1.5/L"
    return prefix + token;
}

SpatialExponentValue GrunObject::interpretRelationship(GrunItem &item)
{
	// 0. data allocation/acquisition
	std::string 				relationship		= item._relationship;			// copy of the item's relationship string (so that the end user's whitespace doesn't get affected)
	SpatialExponentValue		totalRelationshipSV	= SpatialExponentValue::None;	// return value
	std::string 				baseExpr			= "";							// holds the 'base expression' (used to calculate the Spatial 'Rel. Qty')
	std::string					modifierExpr		= "";							// holds the 'modifier expression' (used to compound the Spatial Calculation directly into an Item Qty value - usually this will be for LCIs explicitly using the @ operator)
	bool						isFirstSegment		= true;
	
	// remove whitespace from relationship for processing - whitespace is meaningless, but we want to preserve it in the UI for the end-user
	relationship.erase(std::remove_if(relationship.begin(), relationship.end(), ::isspace), relationship.end());

	// 1. split relationship into implied terms. implied terms are defined by the following regex
	// 		regex looks for: any or no explicit compounding operators, followed by
	//						 any or no numbers, followed by
	//						 any or no GrunObject Tokens
	// dev-note: explicit operators take precedence and if found, define the left-most side of a term
	std::regex splitOnObjectTokens(R"([\/*+\-@]*[0-9.]*[A-Z]*)");
	auto segments_begin		= std::sregex_iterator(relationship.begin(), relationship.end(), splitOnObjectTokens);
	auto segments_end		= std::sregex_iterator();

	// loop through each term
	for (std::sregex_iterator i = segments_begin; i != segments_end; ++i)
	{
		// data acqusition
		std::string	rawSegment			= i->str();								// rawSegment is the current term being processed
		std::string	processedSegment	= injectImplicitOperators(rawSegment);	// processedSegment is the rawSegment (current term being processed) with implicit operators injected into it
		std::string	connection			= "";									// connection is used to insert implicit addition operators between terms in the relationship string

		if (!rawSegment.empty())
		{
			// if this is not the first segment and is a basic term (it starts with a number or uppercase letter), set the connection operator to +
			if (!isFirstSegment)
			{
				if (std::isdigit(rawSegment[0]) || std::isupper(rawSegment[0]))
					connection = "+";
			}

			// check for explicit operators and override the simple + connection operator if needed
			if (rawSegment.starts_with('@'))
			{
				// the @ operator is a compunding explicit operator. It is used to calculate Lineally Centred Items. It allows the end user to directly calculate the Item Qty without using a GrunItem::_itemQuantityFormula. This means it creates a 'modifierExpression' which will NOT be used in the Spatial Value calculation
				std::string val	= rawSegment.substr(1);
				modifierExpr = "/" + val + "+1";
			}
			else if (rawSegment.starts_with('/') || rawSegment.starts_with('*') || rawSegment.starts_with('-') || rawSegment.starts_with('+'))
			{
				// if we're in here it means the user has used an explicit, standard math, operator. When this happens, we have to use the explicit operator OUTSIDE of it's trailing term's parentheses. So the operator must be placed at the left-most character of the term, THEN the opening parenthesis, then the term, then the closing parenthesis
				connection = rawSegment.substr(0,1);	// make the leading, explicit operator the connection character
				processedSegment.erase(0,1);					// pop the leading, explicit operator off the term
			}

			if (!rawSegment.starts_with('@'))
			{
				// join the current segment to the existing baseExpr with the connection operator
				if (!processedSegment.empty())
					baseExpr += connection + "(" + processedSegment + ")";
				isFirstSegment = false;
			}

			// calculate spatial value for this segment
			int segmentExponentTotal = 0;
			std::regex tokenFinder("[LWDAVC]");
			auto t_begin	= std::sregex_iterator(processedSegment.begin(), processedSegment.end(), tokenFinder);
			auto t_end		= std::sregex_iterator();
			for (std::sregex_iterator it = t_begin; it != t_end; ++it)
			{
				// inside a segment, multiplication is the implied operator, so we add the exponent value to the segment's exponent total
				segmentExponentTotal += static_cast<int>(getTokenExponent(it->str()));
			}

			// check for explicit operator - look for it in rawSegment because processedSegment has injected implied operators and may be altered by this point
			// if the user explicitly used an operator, treat the coefficient as a dimension in 3D space, which adds or removes the SpatialExponentValue of the coefficient to the segmentExponentTotal
			size_t opPos = rawSegment.find_first_of("*/+-@");
			bool hasExplicitOperator = (opPos != std::string::npos);
			if (hasExplicitOperator)
			{
				char foundOp = rawSegment[opPos];
				// if (foundOp == '/' || foundOp == '@')
				// {
				// 	// dividing reduces the segmentExponentTotal
				// 	segmentExponentTotal -= 1;
				// 	// std::println("  {:<15} > Reductive operator '{}' detected. Decreasing Spatial Value.", relationship,foundOp);
				// }	
				// else if (foundOp == '*') 
				// {
				// 	// multiplying increases the segmentExponentTotal
				// 	segmentExponentTotal += 1;
				// 	// std::println("  {:<15} > Multiplicative operator '*' detected. Increasing Spatial Value.",relationship);
				// }
				// else 
				// {
				// 	// adding/subtracting do nothing
				// 	// std::println("  {:<15} > Additive operator '{}' detected. Spatial Value remains unchanged.", relationship, foundOp);
				// }
			}

			// the return value takes the highest segmentExponentTotal found from all the segments, ensuring the result is within 0-3
			int safeSegmentExponentTotal = std::clamp(segmentExponentTotal, 0, 3);
			if (safeSegmentExponentTotal > static_cast<int>(totalRelationshipSV)) 
			{
				totalRelationshipSV = static_cast<SpatialExponentValue>(safeSegmentExponentTotal);
			}
		}
	}

	// assign the GrunItem's member values and then return
	item._baseExpression = baseExpr;
	if (!modifierExpr.empty())
	{
		item._interprettedRelationship = "(" + baseExpr + ")" + modifierExpr;
	}
	else
	{
		item._interprettedRelationship = baseExpr;
	}
	item._spatialUnit = totalRelationshipSV;
	return totalRelationshipSV;
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

SpatialExponentValue GrunObject::calculateRelationshipSpatialExponent(const std::string& relationship) const
{
    // ... (This function remains as implemented in the last successful response, handling '*' and '+/-' rules)
    std::regex termSeparatorRegex("([^\\+\\-]+)");
    std::regex dimensionTokenRegex("([LWDAVC])");
    int maxExponent = 0; 

    std::sregex_iterator termIterator(relationship.begin(), relationship.end(), termSeparatorRegex);
    std::sregex_iterator termEnd;

    while (termIterator != termEnd) 
    {
        std::smatch currentTermMatch = *termIterator;
        std::string currentTerm = currentTermMatch.str();
        currentTerm.erase(0, currentTerm.find_first_not_of(" \t\n\r\f\v"));
        currentTerm.erase(currentTerm.find_last_not_of(" \t\n\r\f\v") + 1);

        int currentTermExponent = 0;

        if (currentTerm.find('*') != std::string::npos) 
        {
            // Case A: Multiplication (SUM exponents)
            std::sregex_iterator tokenIterator(currentTerm.begin(), currentTerm.end(), dimensionTokenRegex);
            std::sregex_iterator tokenEnd;

            while (tokenIterator != tokenEnd)
            {
                std::smatch tokenMatch = *tokenIterator;
                std::string token = tokenMatch[1].str(); 
                auto it = propertyDimensions.find(token);
                if (it != propertyDimensions.end()) 
                {
                    currentTermExponent += static_cast<int>(it->second);
                }
                tokenIterator++;
            }
            if (currentTermExponent > static_cast<int>(SpatialExponentValue::Volume)) 
            {
                currentTermExponent = static_cast<int>(SpatialExponentValue::Volume);
            }
        }
        else 
        {
            // Case B: Implicit Addition / Single Dimension (MAX exponent)
            std::sregex_iterator tokenIterator(currentTerm.begin(), currentTerm.end(), dimensionTokenRegex);
            std::sregex_iterator tokenEnd;
            int maxTokenExponent = 0;

            while (tokenIterator != tokenEnd)
            {
                std::smatch tokenMatch = *tokenIterator;
                std::string token = tokenMatch[1].str(); 
                auto it = propertyDimensions.find(token);
                if (it != propertyDimensions.end()) 
                {
                    int tokenExponent = static_cast<int>(it->second);
                    if (tokenExponent > maxTokenExponent)
                    {
                        maxTokenExponent = tokenExponent;
                    }
                }
                tokenIterator++;
            }
            currentTermExponent = maxTokenExponent;
        }
        
        if (currentTermExponent > maxExponent)
        {
            maxExponent = currentTermExponent;
        }
        termIterator++; 
    }
    
    return static_cast<SpatialExponentValue>(maxExponent);
}

/**
 * @brief This function is responsible for all operations involved in determining all values associated with a GrunItem's Spatial Values. It will look at the item's _relationship value and determine all necessary Spatial related values and GrunItem member values when this function is called.
 * @param (GrunItem) - the Grun Item instance to interpret/calculate the Spatial Values for
 * @return (bool) - false if unsuccessful (this will occur if the Grun Item's relationship string is faulty), true if successful.
 */
bool GrunObject::interpretGrunItemSpatialValues(GrunItem &item)
{
	// zero-check and return out false to indicate interpretation was failure
	std::string	relStr			= item._relationship;
	std::string	baseExpr		= "";
	// get the base expression (anything before the last occurence of an @ char, or the whole string)
	auto		atPos			= relStr.find_last_of('@');
	baseExpr = relStr.substr(0, atPos);
	std::string	resultPattern	= "";
	// use regex_replace to find significant characters in the baseExpr and put them in a string														// the replacement pattern used for a regex_replace call
	
	std::string	saneBaseExpr	= std::regex_replace(baseExpr, REGEX_GI_BASEEXPR_SIG_TOKENS_AND_OPS, resultPattern);
	if (saneBaseExpr.empty())
		return false;

	item._baseExpression = baseExpr;
	// data acquisition - make a copy of _relationship because we need to modify it, but preserve the original value
	auto		current			= std::source_location::current();							// for debugging output if needed
	int			spatialAnchor	= 0;

	// we can assume the saneBaseExpr has *something* in it at this point
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
	item._baseExpressionIntprForSU = saneBaseExpr;

	// now we have to calculate what the saneBaseExpr equals in Spatial Unit Value
	// set the saneBaseExpr's total Spatial Value to 0
	int exprTotalSV	= 0;
	
	std::string numericExpr = saneBaseExpr;

	// dev-note: we are assuming saneBaseExpr is a SANITIZED string. If you get unexpected behaviour, you should probably check the value of saneBaseExpr
	// convert all GrunObject Tokens in saneBaseExpr to their Spatial Values
	for (char& c : numericExpr) 
	{
		// convert * operators to +, else convert GrunObject Token to its SpatialExponentValue and cast the number back as a char in the numericExpr string
		if (c == '*')
			c = '+';
		else
			c = static_cast<char>(asInt(getTokenExponent(c)) + 48);	// dev-note: + 48 to correctly cast the returned int value BACK to a char
	}

	item._baseExpressionIntprNumeric = numericExpr;

	// get the largest number in the numericExpr, clamp to min 0 and max 3
	auto digits	= numericExpr 
				| std::views::filter(::isdigit)
            	| std::views::transform([](char c) { return c - '0'; });

	if (!digits.empty())
		spatialAnchor 	= std::clamp(std::ranges::max(digits),0,3);
	item._spatialAnchor	= static_cast<SpatialExponentValue>(spatialAnchor);						  

	// process the operators - loop through the numericExpr char by char with the index value avaiable
	std::string	numericExprResult;
	for (auto [i, c] : std::views::enumerate(numericExpr))
	{
		std::optional<int> r;
		if (c == '+' && i > 0 && i < numericExpr.length())
		{
			const int	lhs	= numericExpr[i-1] - '0';						// dev-note: we must use - '0' here to convert the CHAR 1 into type int 1
			const int	rhs	= numericExpr[i+1] - '0';
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
	
	item._spatialUnit = static_cast<SpatialExponentValue>(spatialUnit);

	// to calculate the Spatial Quantity we must evaluate the GrunItem's base expression
	item._spatialQuantityFormula = convertSpatialQuantitySHNToPEDMAS(item._baseExpression);
	item._spatialQuantityFormula = substituteRelationshipTokens(item._spatialQuantityFormula);
	item._spatialQuantity = evaluateArithmetic(item._spatialQuantityFormula);

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

	// return true to indicate interpretation was a success
	return true;
}

/**
 * @brief Parses a GrunItem's SHN to (hopefully) have correctly interpretted operators and parentheses for PEDMAS arithmetic evaluation
 * @param (std::string) - the shorthand notation as a string to parse and convert
 * @return (std::string) - the converted result string
 */
std::string GrunObject::convertSpatialQuantitySHNToPEDMAS(const std::string &shn)
{
	// zero-check
	if (shn.empty()) return std::string();

	// use regexes to process the shn string in our custom order of precedence to build the formula PEDMAS how we need it
	// dev-note: regex patterns used here are defined in constants at hte top of this file.

	std::string numericForm = shn;
	std::erase_if(numericForm, [](char c) { return std::isspace(static_cast<unsigned char>(c)); });			// strip whitespace
	numericForm	= std::regex_replace(numericForm, REGEX_SHN_TO_PEDMAS_2_NUM_FACTOR_AND_GO_TOKEN, "($1*$2)");
	numericForm = std::regex_replace(numericForm, REGEX_SHN_TO_PEDMAS_4_MISSING_NUMERIC_FACTOR, "(1*");
	numericForm = std::regex_replace(numericForm, REGEX_SHN_TO_PEDMAS_5_IMPLICIT_ADD_OPERATORS, ")+(");

	std::smatch match;
	if (!numericForm.contains('+'))
	{
		// if the numericForm doesn't contain a + operator, the Spatial Value will equal whatever the solitary GrunObject Token represents numerically
		std::regex_search(numericForm, match, REGEX_GO_ALL_TOKENS);
		numericForm = match[0].str();
	}
	else if (std::regex_search(numericForm, match, REGEX_SPATIAL_QTY_SIMPLIFY))
	{
		numericForm = match[0].str();
	}
	return numericForm;
}

bool GrunObject::interpretGrunItemItemQuantity(GrunItem &item)
{
	// zero-check
	if (item._relationship.empty()) return false;

	std::string	parsedRel	= item._relationship;
	std::erase_if(parsedRel, [](char c) { return std::isspace(static_cast<unsigned char>(c)); });			// strip whitespace
	//parsedRel	= std::regex_replace(parsedRel, REGEX_SHN_TO_PEDMAS_0_WRAP_ALL_IN_PARENTHESES, "($1)");
	parsedRel	= std::regex_replace(parsedRel, REGEX_SHN_TO_PEDMAS_1_EXPLICIT_COMBINE_OPERATOR, ")$1(");
	parsedRel	= std::regex_replace(parsedRel, REGEX_SHN_TO_PEDMAS_2_NUM_FACTOR_AND_GO_TOKEN, "($1*$2)");
	parsedRel	= std::regex_replace(parsedRel, REGEX_SHN_TO_PEDMAS_3_AT_OPERATOR, "/$2+1");
	parsedRel	= std::regex_replace(parsedRel, REGEX_SHN_TO_PEDMAS_4_MISSING_NUMERIC_FACTOR, "(1*");
	parsedRel	= std::regex_replace(parsedRel, REGEX_SHN_TO_PEDMAS_5_IMPLICIT_ADD_OPERATORS, ")+(");
	parsedRel	= std::regex_replace(parsedRel, REGEX_SHN_TO_PEDMAS_6_PRIORITIZE_COMBINING_TERMS, "($1)");

	item._interprettedRelationship = substituteRelationshipTokens(parsedRel);
	double	itemQty	= evaluateArithmetic(item._interprettedRelationship);
	// itemQty at this point is only preliminary. check if the item has a value for its _itemQuantityFormula and apply it if it does
	if (!item._itemQuantityFormula.empty())
	{
		std::string	tempForm	= std::to_string(itemQty) + item._itemQuantityFormula;
		itemQty					= evaluateArithmetic(tempForm);
	}
	
	item._itemQuantity = itemQty;
	// std::println("parsedRel: {}",parsedRel);
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
// overload to use char argument
SpatialExponentValue GrunObject::getTokenExponent(char token)
{
	std::string s(1,token);
	return getTokenExponent(s);
}
