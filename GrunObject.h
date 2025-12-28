#include <unordered_map>
#include <string>
#include <vector>
#include "GrunItem.h"

#ifndef CLASS_NAME
#define CLASS_NAME "GrunObject"
#endif

/* GrunObject Class
	dev-todo: 
		1. flesh out GrunItem members (they are the columns from the inventory spreadsheet (Items))
*/

enum class ShapeType
{
	Unknown,
	Rectangle,
	Triangle,
	Circle
};

enum class AreaType
{
	Unknown, 
	Horizontal, 
	Vertical
};

struct TotalAndUnit
{
	double		_total 	= 0.0;
	std::string	_unit	= "unit(s)";
};


// --- CUSTOM FORMATTER FOR TotalAndUnit ---
// This allows TotalAndUnit to be used directly in std::format, std::println, etc.
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
 * @brief Holds 'totals' values about the GrunItems in this GrunObject
 * @property _labourTotal				(std::unordered_map<std::string, TotalsAndUnits>)	first() holds string "Labour", second() holds TotalsAndUnits object with total value and units of measure for the GrunItem
 * @property _materialTotalsSpatial		(std::unordered_map<std::string, TotalsAndUnits>)	first() holds unique GrunItem names found in the GrunObject, second() holds TotalsAndUnits object with total value and units of measure for the GrunItem
 * @property _materialTotalsItemUnit	(std::unordered_map<std::string, TotalsAndUnits>)	first() holds unique GrunItem names found in the GrunObject, second() holds TotalsAndUnits object with total value and units of measure for the GrunItem
 */
struct GrunObjectTotals
{
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
 * @class GrunObject
 * @brief A GrunObject creates a 3D geometric shape object in data (memory).
 * * GrunObject allows developers to create simple, in-memory objects that are based upon 2D geometric shapes (like: Rectangle, Triangle, Circle etc.)
 * * **Key Features**
 * - **Minimal Input** Requires only sizes of the shape's dimensions (x,y,z) and the ShapeType.
 * - **Automatic Attribute Calculations** Automatically calculates additional attributes about the shape (area, volume, aspect ratio, circumference etc.) depending on the ShapeType
*/

class GrunObject
{
	public:
    explicit		GrunObject(const std::string &typeName,
					   const std::string &name,
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
    static int 					asInt(SpatialExponentValue unit);


	bool			addGrunItem(std::string name,					std::string relationship, 
								std::string quantityFormula = "",	std::string units = "unit(s)", 
								std::string primaryLabourFormula = ""
							   );


	int				calculateGrunObjectTotals();
    double 			getAspectRatio();
	std::string		getGrunItemListInfoAsString(const std::string dateFormat = "%d/%m/%Y");
	std::string		getGrunObjectTotalsInfoAsString() const;
	std::string		getObjectName();
	double			getObjectProperty(const std::string propertyName);
	size_t 			removeGrunItem(const std::string& itemName, bool removeAll = false);
	bool			removeGrunItem(size_t index);
	size_t			removeGrunItem(std::vector<size_t> indices);

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
	std::unordered_map<
		std::string, 
		TotalAndUnit> 
			GrunObjectTotals::*					m_totalsPtrs[];

	/**
     * @brief Static constant map that links GrunObject properties (as strings) to their dimensional exponent value.
     * This is the core of the Dimensional Inference Engine (RIE).
     */
    static const	std::unordered_map<std::string, SpatialExponentValue>	propertyDimensions;

	double					applyFormula(double lhs, const std::string &formula, const std::string &itemName, const std::string &type);
	bool					calculateGrunItemData(GrunItem &item);
	SpatialExponentValue	calculateRelationshipSpatialExponent(const std::string& relationship) const;	// class method that determines Spatial Value an item's relationship value is creating
	SpatialExponentValue	calculateGrunItemBaseExpression(GrunItem &item);
	double 					evaluateArithmetic(std::string expression);
	SpatialExponentValue	getTokenExponent(std::string_view token);
	std::string 			substituteRelationshipTokens(const std::string& relationship) const;
	SpatialExponentValue	mapUnitToSpatialExponent(const std::string& unit) const;
	std::string				injectImplicitOperators(std::string &segment);
	SpatialExponentValue	interpretRelationship(GrunItem &item);
};