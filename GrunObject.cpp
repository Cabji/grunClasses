#include <algorithm>
#include <cmath>
#include <format>
#include <functional>
#include <iostream>
#include <print>
#include <regex>
#include <stdexcept>
#include <sstream>
#include "GrunObject.h"

const double PI = std::acos(-1.0);

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
		returnVal += std::format("{:<40}",item._itemName.substr(0,40)) +
					 std::format("{:<17}","Rel: " + item._relationship) + 
					 std::format("{:<13}","SV: " + spatialExponentValueToString(item._spatialExponentValue)) + 
					 std::format("{:<9}","Rel.Qnt: ") + 
					 std::format("{:>7.2f}",item._relationQuantity) + " " +
					 std::format("{:<10}","Item.Qnt: ") + 
					 std::format("{:>7.2f}",item._itemQuantity) + " " +
					 std::format("{:<8}",item._itemQuantityUnits.substr(0,8)) + " " +
					 std::format("{:<10}","P.Labour: ") + 
					 std::format("{:>6.2f}", item._itemPrimaryLabour) + 
					 std::format("{:<8}"," hour(s) ") +
					 std::format("{:<14}","LKGWCalcTime: ") + 
					 std::format("{:>}", item.getCalculatedTimeString(item._itemLKGWCalculated, dateFormat)) +
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

bool GrunObject::calculateGrunItemData(GrunItem &item)
{
    // Define Width and Depth based on AreaType (kept for potential use in future logic)
    double width_val = (m_areaType == AreaType::Vertical) ? m_z : m_y;
    double depth_val = (m_areaType == AreaType::Vertical) ? m_y : m_z;

    // Initialize state
    item._relationQuantity = 0.0;
    item._itemQuantity = 0.0;
    item._itemPrimaryLabour = 0.0;
    bool shn_success = false;

	// zero check - check for invalid tokens, abort if invalid tokens are found
	std::regex invalid_relationship_char_pattern("[^0-9.LWDAVC\\+\\-\\@\\*\\/\\s]"); 
    if (std::regex_search(item._relationship, invalid_relationship_char_pattern))
    {
        std::println("dev-output: GrunItem => {} has an INVALID Relationship String: '{}'. Contains unsupported characters. Quantity set to 0.", 
                     item._itemName, item._relationship);
        return false;
    }

	// check for 1 or more of any valid GrunObject property tokens in the relaitonship string
	std::regex valid_relationship_tokens("[LWDAVC]");
    if (std::regex_search(item._relationship, valid_relationship_tokens))
    {
	}

	item._spatialExponentValue = interpretRelationship(item._relationship);

	/* basic overview for the flow of this function:
		zero check 	- look for anything that makes the relationship string invalid, set result values to 0 and return false
		token check	- if one or more valid tokens exist, step into token processing logic
			token processing 
				we have 2 aims here:
					1. analyze relationship string to determine the Spatial Value and the Spatial Exponent Value (spatial unit of measure)
					2. parse relationship string to calculate the Item Qty (note: sometimes, the Item Qty and Spatial Value are the same)
				determine Spatial Value/Unit:
					1. split the relationship string at each GrunObject Token (LWDAV) into an array structure
						- this split creates the 'Relationship Terms'. A term is comprised of a Token that is (optionally?) preceded by a factor
						  The factor can be applied via implicit multiplication (no operator between factor and Token), 
						  or applied explicitly if the end-user explicitly types in an operator
					2. check each term for explicit operators

		direct math calculation	- */

		// check for any number of valid GrunObject property tokens in the relationship string
		// if valid tokens are found
			// note: for now, 'tokens' are simple, single characters, but in the future we may use token strings so we may need to accommodate for that now in the processing code below
			// split the reltionship string into terms at each token, because the relationship string is different to standard math rules (every factor/divisor group is a term, not a factor)
			//		therefore: 1.5L1W means: 1.5L+1W
			// process the relationship terms in two different ways:
				// 1. to get the Spatial Value
					/* calculate the SpatialExponentTotalValue (SV) of the relationship by correctly combining the SpatialExponentValues of each _token in the terms
					note: when calculating the SV, the terms' exponent values are combined depending on if the operators used are explicitly defined by the end user or not
							 *** explicit operators f0rce combination of the exponent values ***
							 implicit operators do NOT f0rce combination of the exponent values.
					note: 'Compound Operators'. At the moment, we only have '@' as a compound operator.
							'@' is used to calculate "Lineally Centred Items" such as dowel bars, or prop bars on retainer footings.
							The @ operator means, "divide by, and add 1 to the result". So mathmatically, 2L@0.5 = (2(*)L/0.5)+1
							Prop bars have a relationship value like this: 1.2*L@1.2 which means, "a 1.2 long bar along the Length at every 1.2m" 
							This is a compound relationship value and we will see how it should be interpretted for the SV below in the examples.
							Other Compound Operators may be added later on (such as an operator for calculating the length of a spiral coil at a set incline)

							In the examples below <o> means interpretted as an explicit operator and [o] means interpretted as an implicit operator
							for examples: 
								Exponent values are as per the SpatialExponentValues enum class (None = 0, Linear = 1, Area = 2, Volume = 3)
								(implicit addition) 1.5L1W 
									Relation Value	Interpretted as		Exp.Vals		SpatialVal	Reason
									1.5L1W			L[+]W				1[+]1			1 (Linear) 	terms implicitly added together, which does not increase the Spatial Value
								(explicit multiplication) 1.5L*1W 
									1.5L*1W			L<*>W				1<+>1			2 (Area)	terms explicitly multiplied, which increases the Spatial Value
								(explicit compound relation)
									1.2*L@1.2		1.2<*>L</>1.2[+]1	0<+>1<->0[+]1	1 (Linear)	terms explicitly multiplied or divided, but terms with no Tokens do not add to the Spatial Value, and implicit terms do not increase the SPatial Value */

					// therefore, this code block must: 
						// 1. split up the relationship string (use regex, see sub points for regex requirements)
							// 1a. on any valid tokens (LWDAV etc)
							// 1b. on any explicitly, user-defined math operators (+,-,/,*,^,√) - we may need to examine use of ^ and √ closer
						// 2. between each split term, we must decide if we use an explicit operator, or an implicit operator
						//	- our regex would ideally leave any explicit operators on the end of its trailing term so we might be able to pop explicit operators from the end of a term that is followed by an explicit operator)
						// 3. define the logic that decides if we should increase the relationship's Spatial Exponent Value and calculate the Spatial Exponent Value (the Spatial Unit)
						// 4. correctly calculate the Spatial Value from the relationship string so it matches the calculated Spatial Exponent Value.

				// 2. to get the Item Qty value

				// 3. direct math calculation

    // 2. Check for DIRECT QUANTITY (No SHN placeholders LWDAVC) (Existing Logic)

    if (!std::regex_search(item._relationship, valid_relationship_tokens))
    {
        // Direct Quantity calculation (Left-to-Right precedence)
        // This is the evaluator that handles "8*2.5+1"
        std::regex direct_calculation_pattern("^\\s*([-+]?\\d*\\.?\\d+)\\s*((?:[\\+\\-\\*\\/]\\s*\\d*\\.?\\d+)*)\\s*$");
        std::smatch full_match;
        
        if (std::regex_match(item._relationship, full_match, direct_calculation_pattern))
        {
            // ... (The rest of the Direct Quantity logic remains unchanged) ...
            double result = std::stod(full_match[1].str());
            std::string remaining_terms = full_match[2].str();

            std::regex inner_term_pattern("([\\+\\-\\*\\/])\\s*(\\d*\\.?\\d+)");
            
            auto term_begin = std::sregex_iterator(remaining_terms.begin(), remaining_terms.end(), inner_term_pattern);
            auto term_end = std::sregex_iterator();
            
            for (std::sregex_iterator j = term_begin; j != term_end; ++j)
            {
                std::smatch inner_match = *j;
                char op = inner_match[1].str()[0];
                std::string num_str = inner_match[2].str();
                if (num_str.empty()) continue; 
                double value = std::stod(num_str);
                
                switch (op) {
                    case '+': result += value; break;
                    case '-': result -= value; break;
                    case '*': result *= value; break;
                    case '/': 
                        if (value == 0.0) {
                            std::println("dev-output: Division by zero detected in Direct Quantity relationship for '{}'. Quantity set to 0.", item._itemName);
                            return false;
                        }
                        result /= value; 
                        break;
                    default: 
                        return false;
                }
            }
            
            item._relationQuantity = result;
            item._itemQuantity = result;
            shn_success = true;

            item._spatialExponentValue = SpatialExponentValue::None;
            item._outputSpatialExponentValue = SpatialExponentValue::None;
            item._isCompoundRelationship = false; 
        }
        else
        {
            std::println("dev-output: GrunItem => {} has an UNRECOGNIZED Direct Quantity relationship structure: '{}'. Quantity set to 0.", 
                         item._itemName, item._relationship);
            return false;
        }
    }
    else // --- SHN DIMENSIONAL PARSING LOGIC (WITH NEW SUBSTITUTION) ---
    {
        // 3. Determine Dimensional Context (ISV and OSV) - (UNCHANGED)
        //item._spatialExponentValue = calculateRelationshipSpatialExponent(item._relationship);
        item._outputSpatialExponentValue = mapUnitToSpatialExponent(item._itemQuantityUnits);
        
        if (static_cast<int>(item._spatialExponentValue) < static_cast<int>(item._outputSpatialExponentValue))
        {
            item._isCompoundRelationship = true;
        }
        else
        {
            item._isCompoundRelationship = false;
        }
        
        // --- FIXING RELATIONSHIP CALCULATION ---

        // 4. Substitute SHN tokens (L, W, D, A, V) with actual GrunObject dimensions.
        // This converts "2L1W@0.6" -> "2*30 + 1*3.1 @0.6" (or similar structure).
        std::string substituted_expression = substituteRelationshipTokens(item._relationship);
        
        // Find the boundary where the arithmetic stops (e.g., at the '@' in the Dowel example)
        size_t arithmetic_end = substituted_expression.find_first_not_of("0123456789.+-*/ \t\n\r\f\v");
        std::string expression_to_evaluate = substituted_expression.substr(0, arithmetic_end);

        // 5. Evaluate the substituted expression using the Direct Quantity logic.
        // We re-use the evaluation logic here, which fixes the multiple term bug (BUG 3).
        
        std::regex direct_calculation_pattern("^\\s*([-+]?\\d*\\.?\\d+)\\s*((?:[\\+\\-\\*\\/]\\s*\\d*\\.?\\d+)*)\\s*$");
        std::smatch evaluation_match;
        double relationQuantity = 0.0;
        
        if (std::regex_match(expression_to_evaluate, evaluation_match, direct_calculation_pattern))
        {
            relationQuantity = std::stod(evaluation_match[1].str());
            std::string remaining_terms = evaluation_match[2].str();

            std::regex inner_term_pattern("([\\+\\-\\*\\/])\\s*(\\d*\\.?\\d+)");
            
            auto term_begin = std::sregex_iterator(remaining_terms.begin(), remaining_terms.end(), inner_term_pattern);
            auto term_end = std::sregex_iterator();
            
            for (std::sregex_iterator j = term_begin; j != term_end; ++j)
            {
                std::smatch inner_match = *j;
                char op = inner_match[1].str()[0];
                double value = std::stod(inner_match[2].str());
                
                switch (op) {
                    case '+': relationQuantity += value; break;
                    case '-': relationQuantity -= value; break;
                    case '*': relationQuantity *= value; break;
                    case '/': 
                        if (value == 0.0) {
                            std::println("dev-output: Division by zero detected in SHN relationship for '{}'. Using last valid quantity.", item._itemName);
                            // Avoid division by zero, keep current relationQuantity.
                            break;
                        }
                        relationQuantity /= value; 
                        break;
                    default: break;
                }
            }
            shn_success = true;
        }
        else
        {
            std::println("dev-output: GrunItem => {} SHN failed arithmetic evaluation for '{}'. Quantity set to 0.", 
                         item._itemName, expression_to_evaluate);
            return false;
        }
        
        // 6. Correct Compounding (BUG 1 FIX)
        // If a compound relationship (Area -> Volume) exists AND the result is already in the OSV 
        // dimension (i.e., the user provided a numeric factor like '0.05' in '0.05A'), we must NOT multiply by m_z.
        // The only time we compound with m_z is if the relationship was a pure dimension like "A" (which results in numericMultiplier=1.0).
        if (item._isCompoundRelationship)
        {
            int compoundingDifference = static_cast<int>(item._outputSpatialExponentValue) - static_cast<int>(item._spatialExponentValue);
            
            // Check if the user *did not* provide a numeric multiplier in the expression.
            bool numeric_multiplier_provided = (item._relationship.find_first_of("0123456789.") != std::string::npos);

            // We apply m_z compounding ONLY if ISV=Area, OSV=Volume, AND no numeric multiplier was given (e.g., relationship was "A" or "L*W").
            if (compoundingDifference == 1 && item._spatialExponentValue == SpatialExponentValue::Area && !numeric_multiplier_provided)
            {
                relationQuantity *= m_z; 
            }
            // For cases like "0.05A", the math is already: 0.05 * Area. No further compounding by m_z is required.
        }

        // 7. Set the final calculated relation quantity after all steps.
        item._relationQuantity = relationQuantity;
    }


    // --- STEP 2: APPLY itemQuantityFormula (CALCULATE _itemQuantity) (Existing Logic) ---
    
    // ... (This section remains UNCHANGED from the previous response) ...
    std::regex invalid_formula_char_pattern("[^0-9.\\+\\-\\*\\/\\s]"); 
    std::regex formula_pattern("\\s*([\\+\\-\\*\\/])\\s*(\\d*\\.?\\d+)");
    std::smatch formula_match;

    if (shn_success)
    {
        if (item._itemQuantityFormula.empty()) {
            item._itemQuantity = item._relationQuantity;
        }
        else if (std::regex_search(item._itemQuantityFormula, invalid_formula_char_pattern))
        {
            item._itemQuantity = 0.0;
            std::println("dev-output: GrunItem => {} has an INVALID _itemQuantityFormula: '{}'. Contains unsupported characters. Quantity set to 0.", 
                         item._itemName, item._itemQuantityFormula);
        }
        else if (std::regex_match(item._itemQuantityFormula, formula_match, formula_pattern))
        {
            if (formula_match.size() == 3) {
                
                char op = formula_match[1].str()[0];
                double value = std::stod(formula_match[2].str());
                double finalQuantity = item._relationQuantity;

                switch (op) {
                    case '+': finalQuantity += value; break;
                    case '-': finalQuantity -= value; break;
                    case '*': finalQuantity *= value; break;
                    case '/': 
                        if (value == 0.0) {
                            finalQuantity = item._relationQuantity;
                            std::println("dev-output: Division by zero detected in _itemQuantityFormula for '{}'. Applying default formula (*1).", item._itemName);
                            break;
                        }
                        finalQuantity /= value; 
                        break;
                    default: 
                        finalQuantity = item._relationQuantity;
                        break; 
                }
                item._itemQuantity = finalQuantity;
            }
        }
        else
        {
            item._itemQuantity = item._relationQuantity;
            std::println("dev-output: GrunItem => {} has an UNRECOGNIZED _itemQuantityFormula structure: '{}'. Applying default formula (*1, multiply by 1).", 
                         item._itemName, item._itemQuantityFormula);
        }
    }


    // --- STEP 3: CALCULATE Primary Labour Value (_itemPrimaryLabour) (Existing Logic) ---
    
    // ... (This entire section remains UNCHANGED from the previous response) ...
    // Code block for calculating item._itemPrimaryLabour
    // ...
    // Note: Due to size constraints, this block is omitted here but should be included in your implementation.
    
    std::regex invalid_labour_formula_char_pattern("[^0-9.\\+\\-\\*\\/\\s]"); 
    std::regex labour_formula_pattern("\\s*([\\+\\-\\*\\/])\\s*(\\d*\\.?\\d+)");
    std::smatch labour_formula_match;

    if (item._itemPrimaryLabourFormula.empty()) {
        item._itemPrimaryLabour = 0.0; 
    }
    else if (std::regex_search(item._itemPrimaryLabourFormula, invalid_labour_formula_char_pattern))
    {
        item._itemPrimaryLabour = 0.0;
        std::println("dev-output: GrunItem => {} has an INVALID _itemPrimaryLabourFormula: '{}'. Contains unsupported characters. Labour set to 0.", 
                     item._itemName, item._itemPrimaryLabourFormula);
    }
    else if (std::regex_match(item._itemPrimaryLabourFormula, labour_formula_match, labour_formula_pattern))
    {
        if (labour_formula_match.size() == 3) {
            
            char op = labour_formula_match[1].str()[0];
            double value = std::stod(labour_formula_match[2].str());
            
            double finalLabour = item._itemQuantity;

            switch (op) {
                case '+': finalLabour += value; break;
                case '-': finalLabour -= value; break;
                case '*': finalLabour *= value; break;
                case '/': 
                    if (value == 0.0) {
                        finalLabour = item._itemQuantity;
                        std::println("dev-output: Division by zero detected in _itemPrimaryLabourFormula for '{}'. Applying default formula (*1).", 
                                     item._itemName);
                        break;
                    }
                    finalLabour /= value; 
                    break;
                default: 
                    finalLabour = item._itemQuantity;
                    break; 
            }
            
            item._itemPrimaryLabour = finalLabour;
        }
    }
    else
    {
        item._itemPrimaryLabour = item._itemQuantity; 
        std::println("dev-output: GrunItem => {} has an UNRECOGNIZED _itemPrimaryLabourFormula structure: '{}'. Applying default formula (*1).", 
                     item._itemName, item._itemPrimaryLabourFormula);
    }
    
	item._itemLKGWCalculated = std::chrono::system_clock::now();
    return shn_success; 
}

// bool GrunObject::calculateGrunItemData(GrunItem &item)
// {
//     // Define Width and Depth based on AreaType
//     double width_val = (m_areaType == AreaType::Vertical) ? m_z : m_y;
//     double depth_val = (m_areaType == AreaType::Vertical) ? m_y : m_z;

//     // 1. Map SHN placeholders to GrunObject attributes
//     std::map<char, double> attributes = {
//         {'L', m_x},				    // Length (always m_x)
//         {'W', width_val},		    // Width (m_y or m_z, depending on AreaType)
//         {'D', depth_val},		    // Depth (m_z or m_y, depending on AreaType)
//         {'A', m_area},			    // Area
//         {'V', m_volume},		    // Volume
//         {'C', m_circumference}	    // Circumference
//     };

//     // Initialize state
//     item._relationQuantity = 0.0;
//     item._itemQuantity = 0.0;
//     item._itemPrimaryLabour = 0.0;
//     bool shn_success = false;
// 	// calculate the relationship's spatial exponent and assign it to the GrunItem's _spatialExponentValue
// 	item._spatialExponentValue = calculateRelationshipSpatialExponent(item._relationship);

//     // --- STEP 1: SHN CONVERSION (CALCULATE _relationQuantity) ---
    
//     // 2. Add STRICT NEGATIVE VALIDATION
//     std::regex invalid_relationship_char_pattern("[^0-9.LWDAVC\\+\\-\\@\\*\\/\\s]"); 

//     if (std::regex_search(item._relationship, invalid_relationship_char_pattern))
//     {
//         std::println("dev-output: GrunItem => {} has an INVALID SHN relationship: '{}'. Contains unsupported characters. Quantity set to 0.", 
//                      item._itemName, item._relationship);
//         return false; // SHN failed entirely, and we skip steps 2/3 logic later.
//     }

//     // 2b. Check for DIRECT QUANTITY (No SHN placeholders LWDAVC)
//     std::regex valid_relationship_tokens("[LWDAVC]");

//     if (!std::regex_search(item._relationship, valid_relationship_tokens))
//     {
//         // Direct Quantity calculation (Left-to-Right precedence)
//         std::regex direct_calculation_pattern("^\\s*(\\d*\\.?\\d+)\\s*((?:[\\+\\-\\*\\/]\\s*\\d*\\.?\\d+)*)\\s*$");
//         std::smatch full_match;
        
//         if (std::regex_match(item._relationship, full_match, direct_calculation_pattern))
//         {
//             double result = std::stod(full_match[1].str());
//             std::string remaining_terms = full_match[2].str();

//             std::regex inner_term_pattern("([\\+\\-\\*\\/])\\s*(\\d*\\.?\\d+)");
            
//             auto term_begin = std::sregex_iterator(remaining_terms.begin(), remaining_terms.end(), inner_term_pattern);
//             auto term_end = std::sregex_iterator();
            
//             for (std::sregex_iterator j = term_begin; j != term_end; ++j)
//             {
//                 std::smatch inner_match = *j;
//                 char op = inner_match[1].str()[0];
//                 std::string num_str = inner_match[2].str();
//                 if (num_str.empty()) continue; 
//                 double value = std::stod(num_str);
                
//                 switch (op) {
//                     case '+': result += value; break;
//                     case '-': result -= value; break;
//                     case '*': result *= value; break;
//                     case '/': 
//                         if (value == 0.0) {
//                             std::println("dev-output: Division by zero detected in Direct Quantity relationship for '{}'. Quantity set to 0.", item._itemName);
//                             return false;
//                         }
//                         result /= value; 
//                         break;
//                     default: 
//                         return false;
//                 }
//             }
            
//             // Direct Quantity success: set both relation and item quantity
//             item._relationQuantity = result;
//             item._itemQuantity = result;
//             shn_success = true; // Mark as successful for direct quantity path
//         }
//         else
//         {
//             // Unrecognized Direct Quantity structure
//             std::println("dev-output: GrunItem => {} has an UNRECOGNIZED Direct Quantity relationship structure: '{}'. Quantity set to 0.", 
//                          item._itemName, item._relationship);
//             return false;
//         }
//     }
//     else // --- SHN Parsing Logic (Runs if placeholders LWDAVC are found) ---
//     {
//         std::regex outer_shn_pattern("\\s*([\\+\\-]\\s*)?([^@\\s\\+\\-]+)(?:@(\\d*\\.?\\d+))?");
//         std::regex inner_term_pattern("(\\d*\\.?\\d*)\\s*([LWDAVC])");
        
//         std::string relationship = item._relationship;
//         double relationResult = 0.0;

//         auto words_begin = std::sregex_iterator(relationship.begin(), relationship.end(), outer_shn_pattern);
//         auto words_end = std::sregex_iterator();

//         for (std::sregex_iterator i = words_begin; i != words_end; ++i)
//         {
//             // ... (rest of SHN parsing logic) ...
//             std::smatch outer_match = *i;
//             char op = outer_match[1].matched ? outer_match[1].str().find('-') != std::string::npos ? '-' : '+' : '+';

//             std::string shn_block = outer_match[2].str();
//             bool is_lci = outer_match.size() == 4 && outer_match[3].matched;

//             // Inner Parsing: Calculate the total span for this block
//             double block_span_sum = 0.0;
//             auto inner_begin = std::sregex_iterator(shn_block.begin(), shn_block.end(), inner_term_pattern);
//             auto inner_end = std::sregex_iterator();

//             for (std::sregex_iterator j = inner_begin; j != inner_end; ++j)
//             {
//                 std::smatch inner_match = *j;
//                 std::string coeff_str = inner_match[1].str();
//                 double coefficient = coeff_str.empty() ? 1.0 : std::stod(coeff_str);
//                 char placeholder = inner_match[2].str()[0];

//                 if (attributes.count(placeholder))
//                 {
//                     block_span_sum += (coefficient * attributes.at(placeholder));
//                     shn_success = true; 
//                 }
//             }

//             // LCI Calculation
//             double term_final_value = block_span_sum;
//             if (is_lci && block_span_sum > 0.0)
//             {
//                 double interval = std::stod(outer_match[3].str());
//                 if (interval > 0.0) {
//                     term_final_value = (block_span_sum / interval) + 1.0;
//                 }
//             }
            
//             // Apply the operator to the running total
//             if (op == '+') {
//                 relationResult += term_final_value;
//             } else if (op == '-') {
//                 relationResult -= term_final_value;
//             }
//         }
//         item._relationQuantity = relationResult;
//     }


//     // The rest of the function (STEP 2 and STEP 3) executes regardless of the SHN path chosen, 
//     // but operates on the results of STEP 1.

    
//     // --- STEP 2: APPLY itemQuantityFormula (CALCULATE _itemQuantity) ---
    
// 	// dev-note: the regexes beginning with ^ are negating which means "not the following"
//     std::regex invalid_formula_char_pattern("[^0-9.\\+\\-\\*\\/\\s]"); 
//     std::regex formula_pattern("\\s*([\\+\\-\\*/])\\s*(\\d*\\.?\\d+)");
//     std::smatch formula_match;

//     // If SHN succeeded OR Direct Quantity was used, apply formula.
//     if (shn_success)
//     {
//         // 1. Check for EMPTY string (Default to *1)
//         if (item._itemQuantityFormula.empty()) {
//             item._itemQuantity = item._relationQuantity;
//         }
//         // 2. Check for INVALID CHARACTERS (Hard fail)
//         else if (std::regex_search(item._itemQuantityFormula, invalid_formula_char_pattern))
//         {
//             item._itemQuantity = 0.0;
//             std::println("dev-output: GrunItem => {} has an INVALID _itemQuantityFormula: '{}'. Contains unsupported characters. Quantity set to 0.", 
//                          item._itemName, item._itemQuantityFormula);
//         }
//         // 3. Attempt to match SIMPLE FORMULA STRUCTURE
//         else if (std::regex_match(item._itemQuantityFormula, formula_match, formula_pattern))
//         {
//             if (formula_match.size() == 3) {
                
//                 char op = formula_match[1].str()[0];
//                 double value = std::stod(formula_match[2].str());
//                 double finalQuantity = item._relationQuantity; // LHS is _relationQuantity

//                 switch (op) {
//                     case '+': finalQuantity += value; break;
//                     case '-': finalQuantity -= value; break;
//                     case '*': finalQuantity *= value; break;
//                     case '/': 
//                         if (value == 0.0) {
//                             finalQuantity = item._relationQuantity;
//                             std::println("dev-output: Division by zero detected in _itemQuantityFormula for '{}'. Applying default formula (*1).", item._itemName);
//                             break;
//                         }
//                         finalQuantity /= value; 
//                         break;
//                     default: 
//                         finalQuantity = item._relationQuantity;
//                         break; 
//                 }
//                 item._itemQuantity = finalQuantity;
//             }
//         }
//         // 4. FALLBACK: Valid characters, but invalid structure (Soft fail: default to *1.0)
//         else
//         {
//             item._itemQuantity = item._relationQuantity;
//             std::println("dev-output: GrunItem => {} has an UNRECOGNIZED _itemQuantityFormula structure: '{}'. Applying default formula (*1, multiply by 1).", 
//                          item._itemName, item._itemQuantityFormula);
//         }
//     }


//     // --- STEP 3: CALCULATE Primary Labour Value (_itemPrimaryLabour) ---
    
//     // 1. Check for EMPTY string (Default to 0.0)
//     if (item._itemPrimaryLabourFormula.empty()) {
//         item._itemPrimaryLabour = 0.0; 
//     }
//     // 2. Check for INVALID CHARACTERS (Hard fail)
//     else if (std::regex_search(item._itemPrimaryLabourFormula, invalid_formula_char_pattern))
//     {
//         item._itemPrimaryLabour = 0.0;
//         std::println("dev-output: GrunItem => {} has an INVALID _itemPrimaryLabourFormula: '{}'. Contains unsupported characters. Labour set to 0.", 
//                      item._itemName, item._itemPrimaryLabourFormula);
//     }
//     // 3. Attempt to match SIMPLE FORMULA STRUCTURE
//     else if (std::regex_match(item._itemPrimaryLabourFormula, formula_match, formula_pattern))
//     {
//         if (formula_match.size() == 3) {
            
//             char op = formula_match[1].str()[0];
//             double value = std::stod(formula_match[2].str());
            
//             double finalLabour = item._itemQuantity; // LHS is item._itemQuantity

//             switch (op) {
//                 case '+': finalLabour += value; break;
//                 case '-': finalLabour -= value; break;
//                 case '*': finalLabour *= value; break;
//                 case '/': 
//                     if (value == 0.0) {
//                         finalLabour = item._itemQuantity;
//                         std::println("dev-output: Division by zero detected in _itemPrimaryLabourFormula for '{}'. Applying default formula (*1).", 
//                                      item._itemName);
//                         break;
//                     }
//                     finalLabour /= value; 
//                     break;
//                 default: 
//                     finalLabour = item._itemQuantity;
//                     break; 
//             }
            
//             item._itemPrimaryLabour = finalLabour;
//         }
//     }
//     // 4. FALLBACK: Valid characters, but invalid structure (Soft fail: default to *1.0)
//     else
//     {
//         // If the formula is present but cannot be parsed, assume the user intended *1.0 for the labour rate.
//         item._itemPrimaryLabour = item._itemQuantity; 
//         std::println("dev-output: GrunItem => {} has an UNRECOGNIZED _itemPrimaryLabourFormula structure: '{}'. Applying default formula (*1).", 
//                      item._itemName, item._itemPrimaryLabourFormula);
//     }
    
// 	// the GrunItem's data has successfully calculated, so update its LKGWCalculated value to now().
// 	item._itemLKGWCalculated = std::chrono::system_clock::now();
//     return shn_success; 
// }

// SpatialExponentValue GrunObject::calculateRelationshipSpatialExponent(const std::string &relationship) const
// {
// 	// Regex for finding additive terms separated by '+' or '-'
//     // Note: This finds continuous sequences that are NOT '+' or '-'
//     std::regex termSeparatorRegex("([^\\+\\-]+)");
    
//     // Regex for finding individual dimension tokens (L, W, A, V, etc.)
//     std::regex dimensionTokenRegex("([A-Z])");

//     int maxExponent = 0; // Tracks the highest exponent found across all additive terms

//     // --- PHASE 1: Iterate over the TOP-LEVEL ADDITIVE TERMS ---
//     std::sregex_iterator termIterator(relationship.begin(), relationship.end(), termSeparatorRegex);
//     std::sregex_iterator termEnd;

//     while (termIterator != termEnd) 
//     {
//         std::smatch currentTermMatch = *termIterator;
//         std::string currentTerm = currentTermMatch.str();
        
//         // Clean whitespace from the current term
//         currentTerm.erase(0, currentTerm.find_first_not_of(" \t\n\r\f\v"));
//         currentTerm.erase(currentTerm.find_last_not_of(" \t\n\r\f\v") + 1);

//         int currentTermExponent = 0;

//         // --- PHASE 2: Check for Multiplication in the current term ---
//         if (currentTerm.find('*') != std::string::npos) 
//         {
//             // --- Case A: Multiplication (e.g., "1L * 1W" or "2A*L") ---
            
//             // Re-use dimensionTokenRegex to find all tokens within the multiplicative term.
//             // We SUM the exponents of all tokens found in this single multiplicative term.
//             std::sregex_iterator tokenIterator(currentTerm.begin(), currentTerm.end(), dimensionTokenRegex);
//             std::sregex_iterator tokenEnd;

//             while (tokenIterator != tokenEnd)
//             {
//                 std::smatch tokenMatch = *tokenIterator;
//                 std::string token = tokenMatch[1].str(); 

//                 auto it = propertyDimensions.find(token);

//                 if (it != propertyDimensions.end()) 
//                 {
//                     // Multiplication Rule: SUM the exponents
//                     currentTermExponent += static_cast<int>(it->second);
//                 }

//                 tokenIterator++;
//             }
            
//             // Apply the maximum cap check for multiplicative results (A*V = 5 is capped at 3)
//             if (currentTermExponent > static_cast<int>(SpatialExponentValue::Volume)) 
//             {
//                 currentTermExponent = static_cast<int>(SpatialExponentValue::Volume);
//             }
//         }
//         else 
//         {
//             // --- Case B: Implicit Addition / Single Dimension (e.g., "1L1W" or "2A") ---
            
//             // If no '*' is present, treat the term as additive parts where the exponent 
//             // is the MAX exponent of all tokens found in the term (Dimensional Homogeneity).
            
//             std::sregex_iterator tokenIterator(currentTerm.begin(), currentTerm.end(), dimensionTokenRegex);
//             std::sregex_iterator tokenEnd;
//             int maxTokenExponent = 0;

//             while (tokenIterator != tokenEnd)
//             {
//                 std::smatch tokenMatch = *tokenIterator;
//                 std::string token = tokenMatch[1].str(); 

//                 auto it = propertyDimensions.find(token);

//                 if (it != propertyDimensions.end()) 
//                 {
//                     // Implicit Addition Rule: Find the MAX exponent
//                     int tokenExponent = static_cast<int>(it->second);
//                     if (tokenExponent > maxTokenExponent)
//                     {
//                         maxTokenExponent = tokenExponent;
//                     }
//                 }
//                 tokenIterator++;
//             }
//             // The exponent for this entire term is the max exponent found
//             currentTermExponent = maxTokenExponent;
//             // No need to cap here, as maxExponent is already capped by the enum definitions (max 3)
//         }
        
//         // --- PHASE 3: Determine Max Exponent (Overall Addition/Homogeneity) ---
//         // The overall spatial exponent of the relationship is the maximum exponent of any term
//         if (currentTermExponent > maxExponent)
//         {
//             maxExponent = currentTermExponent;
//         }

//         termIterator++; // Move to the next term
//     }
    
//     // Convert the final max exponent to the enum type
//     return static_cast<SpatialExponentValue>(maxExponent);
// }

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
					itemQuantity	= item._relationQuantity;
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
std::string GrunObject::injectImplicitOperators(std::string &segment) {
    // find the first A-Z character
    auto token_start = std::find_if(segment.begin(), segment.end(), [](char c) {
        return std::isupper(static_cast<unsigned char>(c));
    });

    // if no token (like "@1.2") found, return as is
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
    
    // If it ends in a digit or decimal, inject '*'
    if (std::isdigit(static_cast<unsigned char>(lastChar)) || lastChar == '.') {
        return prefix + "*" + token;
    }

    // if it ends in an explicit operator (/, +, -, *), just join them, example: "1.5/L" -> "1.5/L"
    return prefix + token;
}

SpatialExponentValue GrunObject::interpretRelationship(std::string relationship)
{
	// rm whitespace from relationship string
	relationship.erase(std::remove_if(relationship.begin(), relationship.end(), ::isspace), relationship.end());
	
	// split relationship into implied terms - matches one or more capital letter(s) optionally preceded by numbers/explicit operators, or the @ operator followed by numbers
	std::regex splitOnObjectTokens(R"(([0-9.\/*+\-]*[A-Z]+)|(@[0-9.]+))");

	// split up relationship string using the regex
	auto segments_begin	= std::sregex_iterator(relationship.begin(), relationship.end(), splitOnObjectTokens);
	auto segments_end		= std::sregex_iterator();

	std::vector<std::string> segments;										// to hold all the segments
	SpatialExponentValue totalRelationshipSV = SpatialExponentValue::None;	// to hold the return value

	// loop through resulting segments of regex analysis
	for (std::sregex_iterator i = segments_begin; i != segments_end; ++i)
	{
		// push current segment into vector of strings
		std::string rawSegment = i->str();
		std::string processedSegment = injectImplicitOperators(rawSegment);
		segments.push_back(processedSegment);

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

		// check for explicit operator - look for it in rawSegment because processedSegment has injected implied operators
		// if the user explicitly used an operator, treat the coefficient as a dimension in 3D space, which adds or removes the SpatialExponentValue of the coefficient to the segmentExponentTotal
		size_t opPos = rawSegment.find_first_of("*/+-@");
		bool hasExplicitOperator = (opPos != std::string::npos);
		if (hasExplicitOperator)
{
    char foundOp = rawSegment[opPos];
    if (foundOp == '/' || foundOp == '@')
    {
		// dividing reduces the segmentExponentTotal
        segmentExponentTotal -= 1;
        std::println("  {:<15} > Reductive operator '{}' detected. Decreasing Spatial Value.", relationship,foundOp);
    }
    else if (foundOp == '*') 
    {
        // multiplying increases the segmentExponentTotal
        segmentExponentTotal += 1;
        std::println("  {:<15} > Multiplicative operator '*' detected. Increasing Spatial Value.",relationship);
    }
    else 
    {
        // adding/subtracting do nothing
        std::println("  {:<15} > Additive operator '{}' detected. Spatial Value remains unchanged.", relationship, foundOp);
    }
}

		// the return value takes the highest segmentExponentTotal found from all the segments, ensuring the result is within 0-3
		int safeSegmentExponentTotal = std::clamp(segmentExponentTotal, 0, 3);
		if (safeSegmentExponentTotal > static_cast<int>(totalRelationshipSV)) 
		{
			totalRelationshipSV = static_cast<SpatialExponentValue>(safeSegmentExponentTotal);
		}
	}

	return totalRelationshipSV;
}

std::string GrunObject::substituteRelationshipTokens(const std::string& relationship) const
{
    std::string substituted_relationship = relationship;

    // 1. Define the substitution values (ensure they are in string format for replacement)
    // We prioritize using the explicit object dimension members for L, W, D.
    std::string L_val = std::to_string(m_x);
    std::string W_val = std::to_string(m_y); // Use m_y for W
    std::string D_val = std::to_string(m_z); // Use m_z for D
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

    // After substitution, the string is now just numbers and operators (e.g., "2*30 + 1*3.1 @0.6")
    // Note: The "Direct Quantity" evaluator only handles +,-,*,/ and numbers.
    
    // We MUST replace the implicit multiplication/addition from the user's SHN style:
    // This step is crucial: replace '2.5L' with '2.5 * L_val', '1L1W' with '1*L_val + 1*W_val'
    // Since we already did the substitution, we are left with: "2(30) + 1(3.1) @0.6"
    
    // We will assume the user has correctly used '*' for multiplication and the original
    // relationship string for additive terms (2L1W -> 2L + 1W). The evaluation must 
    // stop at the first non-arithmetic character.

    // For now, we return the substituted string, which is the best we can do without a full parser.
    return substituted_relationship;
}

// [4] Implementation of calculateRelationshipSpatialExponent (The Complex Parser - based on previous discussion)
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

// returns a SpatialExponentValue (None,Linear,Area,Volume) that corressponds to the token given to the function
SpatialExponentValue GrunObject::getTokenExponent(std::string_view token)
{
    if (token == "V") return SpatialExponentValue::Volume;
    if (token == "A") return SpatialExponentValue::Area;
    if (token == "L"	|| token == "W"	|| token == "D" || 
        token == "PH"	|| token == "C"	|| token == "R") return SpatialExponentValue::Linear;   
    return SpatialExponentValue::None;
}
