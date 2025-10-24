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
 * @brief Gets a vector of pointers that point to the _members in the GrunObjectTotals instance (m_objectTotals)
 * @return std::vector<std::unordered_map<std::string, TotalAndUnit>*>
 * @note This function simply returns the value from m_objectTotals.getMapPointers, so if you add more members to GrunObjectTotal struct, you have to update GrunObjectTotal::GetMapPointers()
 */
std::vector<std::unordered_map<std::string, TotalAndUnit>*> GrunObject::getTotalPointers()
{
	// get the pointers from the m_objectTotals member
	const auto& returnVal =  m_objectTotals.getMapPointers();
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

std::string GrunObject::getGrunObjectTotals()
{
	// Use std::stringstream to build the output string efficiently
    std::stringstream ss;
    
    // Get the maps and iterate using the index to look up the descriptive name
    std::vector<std::unordered_map<std::string, TotalAndUnit>*> mapPointers = m_objectTotals.getMapPointers();

    ss << std::format("--- GrunObject Totals ({}) ---\n", m_name);
    
    for (size_t i = 0; i < mapPointers.size(); ++i)
    {
        const auto* totalMap = mapPointers[i];
        
        // 1. Output the section header (using the descriptive name)
        ss << std::format("Map [{}]: {}\n", i, m_objectTotals.getMapName(i));

        if (totalMap->empty())
        {
            ss << "\t(No aggregated totals in this section)\n";
        }
        else
        {
            // 2. Loop through the key-value pairs in the unordered_map
            // The value 'TotalAndUnit' will be automatically formatted by the custom formatter
            for (const auto& pair : *totalMap)
            {
                // pair.first is std::string (the item name)
                // pair.second is TotalAndUnit (the value)
                ss << std::format("\t- {:20}: {}\n", pair.first, pair.second);
            }
        }
        ss << "\n";
    }

    return ss.str();
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
bool GrunObject::calculateGrunItemData(GrunItem &item)
{
    // Define Width and Depth based on AreaType
    double width_val = (m_areaType == AreaType::Vertical) ? m_z : m_y;
    double depth_val = (m_areaType == AreaType::Vertical) ? m_y : m_z;

    // 1. Map SHN placeholders to GrunObject attributes
    std::map<char, double> attributes = {
        {'L', m_x},				    // Length (always m_x)
        {'W', width_val},		    // Width (m_y or m_z, depending on AreaType)
        {'D', depth_val},		    // Depth (m_z or m_y, depending on AreaType)
        {'A', m_area},			    // Area
        {'V', m_volume},		    // Volume
        {'C', m_circumference}	    // Circumference
    };

    // Initialize state
    item._relationQuantity = 0.0;
    item._itemQuantity = 0.0;
    item._itemPrimaryLabour = 0.0;
    bool shn_success = false;


    // --- STEP 1: SHN CONVERSION (CALCULATE _relationQuantity) ---
    
    // 2. Add STRICT NEGATIVE VALIDATION
    std::regex invalid_shn_char_pattern("[^0-9.LWDAVC\\+\\-\\@\\*\\/\\s]"); 

    if (std::regex_search(item._relationship, invalid_shn_char_pattern))
    {
        std::println("dev-output: GrunItem => {} has an INVALID SHN relationship: '{}'. Contains unsupported characters. Quantity set to 0.", 
                     item._itemName, item._relationship);
        return false; // SHN failed entirely, and we skip steps 2/3 logic later.
    }

    // 2b. Check for DIRECT QUANTITY (No SHN placeholders LWDAVC)
    std::regex shn_placeholder_pattern("[LWDAVC]");

    if (!std::regex_search(item._relationship, shn_placeholder_pattern))
    {
        // Direct Quantity calculation (Left-to-Right precedence)
        std::regex direct_calculation_pattern("^\\s*(\\d*\\.?\\d+)\\s*((?:[\\+\\-\\*\\/]\\s*\\d*\\.?\\d+)*)\\s*$");
        std::smatch full_match;
        
        if (std::regex_match(item._relationship, full_match, direct_calculation_pattern))
        {
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
            
            // Direct Quantity success: set both relation and item quantity
            item._relationQuantity = result;
            item._itemQuantity = result;
            shn_success = true; // Mark as successful for direct quantity path
        }
        else
        {
            // Unrecognized Direct Quantity structure
            std::println("dev-output: GrunItem => {} has an UNRECOGNIZED Direct Quantity relationship structure: '{}'. Quantity set to 0.", 
                         item._itemName, item._relationship);
            return false;
        }
    }
    else // --- SHN Parsing Logic (Runs if placeholders LWDAVC are found) ---
    {
        std::regex outer_shn_pattern("\\s*([\\+\\-]\\s*)?([^@\\s\\+\\-]+)(?:@(\\d*\\.?\\d+))?");
        std::regex inner_term_pattern("(\\d*\\.?\\d*)\\s*([LWDAVC])");
        
        std::string relationship = item._relationship;
        double relationResult = 0.0;

        auto words_begin = std::sregex_iterator(relationship.begin(), relationship.end(), outer_shn_pattern);
        auto words_end = std::sregex_iterator();

        for (std::sregex_iterator i = words_begin; i != words_end; ++i)
        {
            // ... (rest of SHN parsing logic) ...
            std::smatch outer_match = *i;
            char op = outer_match[1].matched ? outer_match[1].str().find('-') != std::string::npos ? '-' : '+' : '+';

            std::string shn_block = outer_match[2].str();
            bool is_lci = outer_match.size() == 4 && outer_match[3].matched;

            // Inner Parsing: Calculate the total span for this block
            double block_span_sum = 0.0;
            auto inner_begin = std::sregex_iterator(shn_block.begin(), shn_block.end(), inner_term_pattern);
            auto inner_end = std::sregex_iterator();

            for (std::sregex_iterator j = inner_begin; j != inner_end; ++j)
            {
                std::smatch inner_match = *j;
                std::string coeff_str = inner_match[1].str();
                double coefficient = coeff_str.empty() ? 1.0 : std::stod(coeff_str);
                char placeholder = inner_match[2].str()[0];

                if (attributes.count(placeholder))
                {
                    block_span_sum += (coefficient * attributes.at(placeholder));
                    shn_success = true; 
                }
            }

            // LCI Calculation
            double term_final_value = block_span_sum;
            if (is_lci && block_span_sum > 0.0)
            {
                double interval = std::stod(outer_match[3].str());
                if (interval > 0.0) {
                    term_final_value = (block_span_sum / interval) + 1.0;
                }
            }
            
            // Apply the operator to the running total
            if (op == '+') {
                relationResult += term_final_value;
            } else if (op == '-') {
                relationResult -= term_final_value;
            }
        }
        item._relationQuantity = relationResult;
    }


    // The rest of the function (STEP 2 and STEP 3) executes regardless of the SHN path chosen, 
    // but operates on the results of STEP 1.

    
    // --- STEP 2: APPLY itemQuantityFormula (CALCULATE _itemQuantity) ---
    
	// dev-note: the regexes beginning with ^ are negating which means "not the following"
    std::regex invalid_formula_char_pattern("[^0-9.\\+\\-\\*\\/\\s]"); 
    std::regex formula_pattern("\\s*([\\+\\-\\*/])\\s*(\\d*\\.?\\d+)");
    std::smatch formula_match;

    // If SHN succeeded OR Direct Quantity was used, apply formula.
    if (shn_success)
    {
        // 1. Check for EMPTY string (Default to *1)
        if (item._itemQuantityFormula.empty()) {
            item._itemQuantity = item._relationQuantity;
        }
        // 2. Check for INVALID CHARACTERS (Hard fail)
        else if (std::regex_search(item._itemQuantityFormula, invalid_formula_char_pattern))
        {
            item._itemQuantity = 0.0;
            std::println("dev-output: GrunItem => {} has an INVALID _itemQuantityFormula: '{}'. Contains unsupported characters. Quantity set to 0.", 
                         item._itemName, item._itemQuantityFormula);
        }
        // 3. Attempt to match SIMPLE FORMULA STRUCTURE
        else if (std::regex_match(item._itemQuantityFormula, formula_match, formula_pattern))
        {
            if (formula_match.size() == 3) {
                
                char op = formula_match[1].str()[0];
                double value = std::stod(formula_match[2].str());
                double finalQuantity = item._relationQuantity; // LHS is _relationQuantity

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
        // 4. FALLBACK: Valid characters, but invalid structure (Soft fail: default to *1.0)
        else
        {
            item._itemQuantity = item._relationQuantity;
            std::println("dev-output: GrunItem => {} has an UNRECOGNIZED _itemQuantityFormula structure: '{}'. Applying default formula (*1, multiply by 1).", 
                         item._itemName, item._itemQuantityFormula);
        }
    }


    // --- STEP 3: CALCULATE Primary Labour Value (_itemPrimaryLabour) ---
    
    // 1. Check for EMPTY string (Default to 0.0)
    if (item._itemPrimaryLabourFormula.empty()) {
        item._itemPrimaryLabour = 0.0; 
    }
    // 2. Check for INVALID CHARACTERS (Hard fail)
    else if (std::regex_search(item._itemPrimaryLabourFormula, invalid_formula_char_pattern))
    {
        item._itemPrimaryLabour = 0.0;
        std::println("dev-output: GrunItem => {} has an INVALID _itemPrimaryLabourFormula: '{}'. Contains unsupported characters. Labour set to 0.", 
                     item._itemName, item._itemPrimaryLabourFormula);
    }
    // 3. Attempt to match SIMPLE FORMULA STRUCTURE
    else if (std::regex_match(item._itemPrimaryLabourFormula, formula_match, formula_pattern))
    {
        if (formula_match.size() == 3) {
            
            char op = formula_match[1].str()[0];
            double value = std::stod(formula_match[2].str());
            
            double finalLabour = item._itemQuantity; // LHS is item._itemQuantity

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
    // 4. FALLBACK: Valid characters, but invalid structure (Soft fail: default to *1.0)
    else
    {
        // If the formula is present but cannot be parsed, assume the user intended *1.0 for the labour rate.
        item._itemPrimaryLabour = item._itemQuantity; 
        std::println("dev-output: GrunItem => {} has an UNRECOGNIZED _itemPrimaryLabourFormula structure: '{}'. Applying default formula (*1).", 
                     item._itemName, item._itemPrimaryLabourFormula);
    }
    
	// the GrunItem's data has successfully calculated, so update its LKGWCalculated value to now().
	item._itemLKGWCalculated = std::chrono::system_clock::now();
    return shn_success; 
}

/**
 * @brief	Calculates the Totals data for the GrunObject.
 * @note	This function utilises a pointer to look at the members of the m_objectTotals object. There are 3 std::unordered_maps in the TotalsAndUnits struct and in this function ptrToAnUnorderedMap points to them consecutively through a for loop clause.
 * @warning	If the design of TotalsAndUnits struct changes, you must ensure contiguity of these unordered_maps members for this function to work correctly.
 */
int GrunObject::calculateGrunObjectTotals()
{
	if (m_items.empty()) { return 0; }

	// local data
	int		returnVal = 0;
	// create a psuedo-type 'MapPointer' and create a vector of MapPointers, load it with pointers to the members in m_objectTotals
    using MapPointer = std::unordered_map<std::string, TotalAndUnit>*;
    std::vector<MapPointer> ptrsToTotalsMembers = m_objectTotals.getMapPointers();
    
	// define the _source_ of Quantity and Unit for each destination map.
    // The order MUST match the order returned by GrunObjectTotals::getMapPointers().
    // We use C++ Pointer-to-Member syntax (e.g., &GrunItem::_itemQuantity) 
    // to dynamically access the correct member for each map.

	// create a psuedo-type 'AccessorPair' and load it with 
    using AccessorPair = std::pair<double GrunItem::*, std::string GrunItem::*>;
    std::array<AccessorPair, 3> accessors = {
        // Map 0 (_labourTotal) -> Quantity from _primaryLabourQuantity, Unit from _primaryLabourUnits
        std::make_pair(&GrunItem::_itemPrimaryLabour, &GrunItem::_itemPrimaryLabourUnits),
        
        // Map 1 (_materialTotalsSpatial) -> Quantity from _relationQuantity, Unit from _relationship
		// dev-note: this assignment is incorrect for now. i have to implement auto detection of spatial units based on _relationship expression
        std::make_pair(&GrunItem::_relationQuantity, &GrunItem::_relationship),
        
        // Map 2 (_materialTotalsItemUnit) -> Quantity from _itemQuantity, Unit from _itemUnits
        std::make_pair(&GrunItem::_itemQuantity, &GrunItem::_itemQuantityUnits)
    };

	// Safety check: ensure the number of accessors matches the number of maps
    if (accessors.size() != ptrsToTotalsMembers.size())
    {
        // This is a configuration error, should ideally be caught at compile time.
        // For runtime safety, return an error code or throw.
        return 1; // Return non-zero for error
    }

	// loop the m_items vector and for each item, decide if the GrunItem's _itemName needs a new entry in the m_objectTotals' unordered_map	members 
	for (const auto& item : m_items)
	{
		// loop over the destination maps (we use the for loop's index value (i) to choose the correct accessor)
		for (size_t i = 0; i < ptrsToTotalsMembers.size(); i++)
		{
			MapPointer ptrToMap = ptrsToTotalsMembers[i];

			// extract the correct accessor pair
			double 			GrunItem::*	quantityMember	= accessors[i].first;
			std::string		GrunItem::*	unitMember		= accessors[i].second;

			// get the values from the current item
			double						itemQuantity	= item.*quantityMember;
	const	std::string&				itemUnit		= item.*unitMember; 
	const	std::string&				itemName		= item._itemName;

			// skip processing if the quantity is zero
			if (itemQuantity == 0.0) continue;

			// Use operator[] to get a reference to the entry. 
            // If the key (itemName) doesn't exist, it inserts a new TotalAndUnit 
            // with default values (total=0.0, unit="").
            TotalAndUnit& currentEntry = (*ptrToMap)[itemName];
			
			// IF the unit is currently empty, set it from the item.
            // This assumes the first time an item is encountered, its unit is the definitive one.
            if (currentEntry._unit.empty())
            {
                currentEntry._unit = itemUnit;
            }
            // ELSE IF the unit is different, you might want to log a warning or throw an error
            // because you cannot aggregate different units (e.g., "m" and "cm").
            else if (currentEntry._unit != itemUnit && !itemUnit.empty())
            {
				std::println("Warning: Attemping to combine '{}' with unit '{}' into an existing entry with unit'{}'. Skipping unit update.", 
							itemName, itemUnit, currentEntry._unit
				);

            }
			// Always add the quantity to the total
            currentEntry._total += itemQuantity;
		}
	}

	return 0;
}