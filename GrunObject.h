#include <format>
#include <optional>
#include <unordered_map>
#include <string>
#include <vector>
#include "GrunItem.h"

#ifndef CLASS_NAME
#define CLASS_NAME "GrunObject"
#endif

/**
 * @brief Assigns a numeric 'exponent' value to S.H.N. GrunObject Tokens. If you add additional GrunObject properties/tokens you would need to assign them an exponent value in this enum class' definition.
 */
enum class GrunObjectMemberSUV 
{
	L	= 1,
	W	= 1,
	D	= 1,
	A	= 2,
	V	= 3,
	C	= 1,
	R	= 1
};

/**
 * @brief Defines known, geometric shape types that a GrunObject can be assigned as, to decide how to calculate its dimensions. If you add more geometric shape types to the GrunObject class, you need to define them in this enum class definition.
 */
enum class ShapeType
{
	Unknown,
	Rectangle,
	Triangle,
	Circle
};

/**
 * @brief Defines which 'plane' the Area of a GrunObject should be calculated on. Horizontal is like a floor, Vertical is like a wall.
 */
enum class AreaType
{
	Unknown, 
	Horizontal, 
	Vertical
};

/**
 * @brief Holds a "Total" value based on "unit". This is used for totalling calculated quantities for all GrunItems in a GrunObject that have a common unit of measure.
 */
struct TotalAndUnit
{
	double		_total 	= 0.0;
	std::string	_unit	= "unit(s)";
};

struct ItemAndTotal
{
	std::string		itemName	= "";
	TotalAndUnit	itemTotal;
};

/**
 * @brief Holds 'totals' values about the GrunItems in this GrunObject
 * @property _labourTotal				(std::unordered_map<std::string, TotalsAndUnits>)	first() holds string "Labour", second() holds TotalsAndUnits object with total value and units of measure for the GrunItem
 * @property _materialTotalsSpatial		(std::unordered_map<std::string, TotalsAndUnits>)	first() holds unique GrunItem names found in the GrunObject, second() holds TotalsAndUnits object with total value and units of measure for the GrunItem
 * @property _materialTotalsItemUnit	(std::unordered_map<std::string, TotalsAndUnits>)	first() holds unique GrunItem names found in the GrunObject, second() holds TotalsAndUnits object with total value and units of measure for the GrunItem
 */
struct GrunObjectTotals
{
	friend class GrunObject;

	private: 
	// maps for storing aggregated totals data about a GrunObject's GrunItems
	std::unordered_map<std::string, TotalAndUnit>	_labourTotal;
	std::unordered_map<std::string, TotalAndUnit>	_materialTotalsSpatial;
	std::unordered_map<std::string, TotalAndUnit>	_materialTotalsItemUnit;

	public:
	using	GrunObjectTotalsPtrs = std::unordered_map<std::string, TotalAndUnit> GrunObjectTotals::*; 
	static	constexpr	GrunObjectTotalsPtrs		TOTALS_PTRS[] =
	{
		&GrunObjectTotals::_labourTotal,
		&GrunObjectTotals::_materialTotalsSpatial, 
		&GrunObjectTotals::_materialTotalsItemUnit
	};


	/**
     * @brief Returns the number of total maps available.
     * @return The size of the MAP_PTRS array.
     */
    static constexpr size_t getMapCount() {
        return std::size(TOTALS_PTRS);
    }

	/**
     * @brief Returns a descriptive name for the map at a given index.
     * @param index The index of the map (0, 1, or 2).
     * @return The descriptive name.
     */
    std::string getMapName(size_t index) const
    {
        switch(index)
        {
            case 0: return "Labour Total (Primary Labour)";
            case 1: return "Material Totals (Spatial Relation)";
            case 2: return "Material Totals (Item Unit Quantity)";
            default: return "Unknown Map - you may want to check GrunObjectTotals::getMapPointers() and update it!";
        }
    }
};

/**
 * @brief	A struct that holds vetted data about a GrunObject's child items.
 * @note	This is used by the getter function that sends data about the object's children out to calling code.
 */
struct GrunItemSummary
{
	std::string	name;
	double		quantity;
	double		cost;
};

/**
 * @struct RelationshipSearchResult
 * @brief A structure that houses index location results from a relationship/comment search on GrunItem's RelationshipValues vector
 */
struct RelationshipSearchResult
{
	size_t	itemIndex;			// index in GrunObject::m_items
	size_t	relationshipIndex;	// index in GrunItem::_itemCoreValues
};

// set the ShapeType from a string - NOTE: we must ensure the string is sanitized to all lowercase if the user inputs this value
inline ShapeType shapeTypeFromString(const std::string &s)
{
    if (s == "rectangle") return ShapeType::Rectangle;
    if (s == "triangle")  return ShapeType::Triangle;
    if (s == "circle")    return ShapeType::Circle;
    return ShapeType::Unknown;
}

// return a string value for whatever ShapeType the object is
inline std::string shapeTypeToString(ShapeType t)
{
    switch (t)
    {
        case ShapeType::Rectangle: return "rectangle";
        case ShapeType::Triangle:  return "triangle";
        case ShapeType::Circle:    return "circle";
        default:                   return "unknown";
    }
}

// set the AreaType from a string - NOTE: we must ensure the string is sanitized to all lowercase if the user inputs this value
inline AreaType areaTypeFromString(const std::string &s)
{
    if (s == "horizontal") return AreaType::Horizontal;
    if (s == "vertical")  return AreaType::Vertical;
    return AreaType::Unknown;
}

// return a string value for whatever AreaType the object has
inline std::string areaTypeToString(AreaType t)
{
    switch (t)
    {
        case AreaType::Horizontal: return "horzintal";
        case AreaType::Vertical:   return "vertical";
        default:                   return "unknown";
    }
}

/**
 * @brief Custom formatter for TotalAndUnit struct. Allows you to easily output TotalAndUnit objects using std::println().
 */
template <>
struct std::formatter<TotalAndUnit> : std::formatter<std::string> 
{
    // Optionally hold formatting options here if you want to support {:#}, {:.2f}, etc.
    // For simplicity, we delegate parsing to the string formatter.

    auto format(const TotalAndUnit& t, std::format_context& ctx) const 
    {
        // Format the object into a simple string: "Total (Unit)"
        std::string s = std::format("{:.3f} {}", t._total, t._unit); // Using {:.3f} for double precision
        
        // Pass the formatted string to the underlying string formatter
        return std::formatter<std::string>::format(s, ctx);
    }
};

/**
 * @brief Custom formatter for ItemAndTotal struct. Allows you to easily output ItemAndTotal objects using std::println().
 */
template <>
struct std::formatter<ItemAndTotal> : std::formatter<std::string> 
{
    // Optionally hold formatting options here if you want to support {:#}, {:.2f}, etc.
    // For simplicity, we delegate parsing to the string formatter.

    auto format(const ItemAndTotal& t, std::format_context& ctx) const 
    {
        // Format the object into a simple string: "Total (Unit)"
        std::string s = std::format("{}: {}", t.itemName, t.itemTotal); // itemTotal output should be handled by the TotalAndUnit custom formatter, yes?
        
        // Pass the formatted string to the underlying string formatter
        return std::formatter<std::string>::format(s, ctx);
    }
};

/**
 * @brief Custom formatter for RelationshipSearchResult struct. Allows you to output a RelationshipSearchResult easily using std::println().
 */
template <>
struct std::formatter<RelationshipSearchResult>
{
	// 1. Tell the compiler we don't need any special flags (like :.2f)
    constexpr auto parse(std::format_parse_context& ctx) {
        return ctx.begin();
    }

    auto format(const RelationshipSearchResult& r, std::format_context& ctx) const 
    {
        // write to the context's output buffer
        return std::format_to(ctx.out(), "[Item: {} => Relationship: {}]", r.itemIndex, r.relationshipIndex);
    }
};

/** 
 * @class GrunObject
 * @brief A GrunObject creates a representation of a 3D geometric shape object (a solid) in data (memory).\n
 * GrunObject allows developers to create simple, in-memory objects that are based upon 2D geometric shapes (like: Rectangle, Triangle, Circle etc.)\n
 * Key Features {key-features}
 * - Minimal Input: Requires only sizes of the shape's dimensions (x,y,z) and the ShapeType.
 * - Automatic Attribute Calculations: Automatically calculates additional attributes about the shape (area, volume, aspect ratio, circumference etc.) depending on the ShapeType
*/
class GrunObject
{
	public:
	GrunObject()	= default;

    explicit		GrunObject(const std::string &shapeType = "unknown",
					   const std::string &name = "DefaultName",
                       double m_x = 0.0,
                       double m_y = 0.0, 
					   double m_z = 0.0,
					   const std::string &areaType = "horizontal",
					   const std::string &stage = "");

    /**
     * @brief Gets the SpatialExponentValue (0, 1, 2, or 3) for a given property name token.
     * @param propertyName The property name token (e.g., "L", "A", "Volume").
     * @return The corresponding SpatialExponentValue enum. Defaults to Unitless (0) if not found.
     */
    static SpatialExponentValue getSpatialUnit(const std::string& propertyName);

    /**
     * @brief Converts a SpatialExponentValue enum to its underlying integer value (0, 1, 2, or 3).
     * @param unit The SpatialExponentValue to convert.
     * @return The integer exponent value.
     */
    static int 			asInt(SpatialExponentValue unit);


	bool							addGrunItem(				std::string 				name,
																std::string 				relationship 			= "", 
																std::string 				relComment 				= "",
																std::string 				quantityFormula 		= "",
																std::string 				units 					= "unit(s)", 
																std::string 				primaryLabourFormula 	= ""
												);

	size_t							addGrunItemRelationship(	GrunItem&					item,
																const std::string&			relationship			= "", 
																const std::string&			relComment				= ""
															);
	size_t							rmGrunItemRelationship(		GrunItem&					item,
																const std::string& 			relationship			= "",
																const std::string& 			relComment				= ""
															);


	std::vector
		<size_t>					findGrunItemByItemName(		const std::string&			findItemName, 
																bool						useExactSearch			= true
														) const;
	std::vector
		<RelationshipSearchResult>	findRelationshipByStrings(	const size_t 				itemIndex,
																const std::string&			relationship,
																const std::string&			relComment,
																const bool 					useExactSearch 			= true
															) const;
	std::vector
		<RelationshipSearchResult>	findRelationshipByStrings(	const std::vector<size_t>	itemIndices,
																const std::string& 			relationship,
																const std::string& 			relComment,
																const bool 					useExactSearch			= true
															) const;
	int								calculateGrunObjectTotals();
    double 							getAspectRatio();
	GrunItem&						getGrunItemByIndex(			const int 					index);
	std::string						getGrunItemListInfoAsString(const std::string 			dateFormat				= "%d/%m/%Y");
	std::string						getGrunObjectTotalsInfoAsString() const;
	std::string						getObjectName();
	ItemAndTotal					getItemQtyTotal(			const size_t& 				index, 
																const bool& 				getRounded				= false
													);
	std::vector
		<ItemAndTotal>				getItemQtyTotals(			const bool& 				getRounded				= false);
	double							getObjectProperty(			const std::string			propertyName);
	size_t							getTotalOfGrunItems();
	size_t 							removeGrunItem(				const std::string&			itemName,
																bool						removeAll				= false);
	bool							removeGrunItem(				size_t						index);
	size_t							removeGrunItem(				std::vector<size_t>			indices);
	size_t							size();
	
	private:
	std::string									m_name;					// the GrunObject's name
	std::string									m_stage;				// the GrunObject's stage assignment
	ShapeType									m_type;					// the basic 2D geometric shape this object is based upon (Rectangle, Triangle etc.)
	double										m_x;					// typically, length
	double										m_y;					// typically, width
	double										m_z;					// typically, depth
	AreaType									m_areaType;				// decides the way an object's 'Area' ought to be calculated (L*W vs L*D)
	double										m_area;
	double										m_volume;
	double										m_aspectRatio;			// aspect ratio of the shape if valid. Uses m_areaType to decide which dimensions to use to obtain the aspect ratio.
	double										m_circumference;		// the circumference of the shape if it's a Circle.
	std::vector<GrunItem>						m_items;				// std::vector of GrunItems associated to the GrunObject
	GrunObjectTotals							m_objectTotals;			// an object that holds Totals data about the GrunItems in this GrunObject
	
	using MapMemberPtr = 	std::unordered_map<std::string, TotalAndUnit>GrunObjectTotals::*;
	static constexpr 
	std::array
		<
		MapMemberPtr, 
		GrunObjectTotals::getMapCount()
	>											m_totalsPtrs	= 	{
																		&GrunObjectTotals::_labourTotal,
																		&GrunObjectTotals::_materialTotalsSpatial,
																		&GrunObjectTotals::_materialTotalsItemUnit
																	};

	/**
     * @brief Static constant map that links GrunObject properties (as strings) to their dimensional exponent value.
     * This is the core of the Dimensional Inference Engine (RIE).
     */
    static const	std::unordered_map<std::string, SpatialExponentValue>	propertyDimensions;

	double					applyFormula(double lhs, const std::string &formula);

	bool					interpretGrunItemSpatialValues(GrunItem &item);
	std::string				convertSpatialQuantitySHNToPEDMAS(const std::string &shn);
	bool					interpretGrunItemItemQuantity(GrunItem &item);
	bool					calculateGrunItemData(GrunItem &item);
	double 					evaluateArithmetic(std::string expression);
	SpatialExponentValue	getTokenExponent(std::string_view token);
	SpatialExponentValue	getTokenExponent(char token);
	std::string 			substituteRelationshipTokens(const std::string& relationship) const;
	SpatialExponentValue	mapUnitToSpatialExponent(const std::string& unit) const;
	std::string 			spatialUnitToString(SpatialExponentValue unit);
	std::string 			spatialUnitToSuffix(SpatialExponentValue unit);
};