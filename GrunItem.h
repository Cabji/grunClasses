#pragma once
#include <chrono>
#include <ctime>
#include <functional>

#ifndef CLASS_NAME
#define CLASS_NAME "GrunItem"
#endif

// used for calculating what spatial unit a GrunItem relationship should result in
enum class SpatialExponentValue {
    None 	= 0,
    Linear  = 1,
    Area    = 2,
    Volume  = 3
};

/**
 * @brief Converts a SpatialExponentValue enum to its human-readable string representation.
 * @param exponent - The enum value to convert.
 * @return std::string - The corresponding string ("None", "Linear", "Area", "Volume").
 */
std::string spatialExponentValueToString(SpatialExponentValue exponent);

/**
 * @brief RelationshipValues struct are the members in Grunitem that are all linked together via the Item's relationship.
 * @note Because we need to support Items having multiple relationships, this struct will help us keep all members that are influenced by the item's relationship in a neat unit.
 */
struct RelationshipValues
{
	std::string				relationship;					// the relationship string
	std::string				relComment;						// the relationship comment string
	// all values below here are interpretted and/or calculated from the relationship string

	std::string				baseExpression;					// the full base expression of the relationship (used for spatial quantity calculations)
	std::string				baseExpressionIntrpForSU;		// the base expression after simplification
	std::string				baseExpressionIntrpNumeric;		// the base expression, after simplification, with tokens switched out for numeric SPatial Exponent values
	SpatialExponentValue	spatialAnchor;					// the Spatial Anchor for this relationship only
	std::string				spatialQuantityFormula;			// ?
	double					spatialQuantity;				// the Spatial Quantity for this relationship only
	SpatialExponentValue	spatialUnit;					// the Spatial Unit for this relationship only
	std::string				interprettedRelationship;		// the full relationship string, after intepretation (with numeric values?)
	double					itemQuantity;					// the item quantity for this relationship only

};

/**
 * @brief GrunItem is an instance of an Item from the inventory (eg: a Material, Service, or some form of Secondary Labour). GrunItem will generally have a 'relationship' to its 'owner' GrunObject, or GrunItem can simply have an incidental amount as its relationship, like: 1 to add 1 roll of tie wire to the job, or 8 to set 8 hours of buffer labour to the job.
 * @note mandatory info about a GrunItem to create a fully working instance of it is: 
	- _itemName, _relationship, _itemQuantityFormula, _itemUnits, _itemPrimaryLabourFormula
	
*	Below is a logical walkthrough of how a GrunItem's data is built from just the SHN and its owner GrunObject
*		_relationship:				ShortHand Notation (SHN) for how the GrunItem relates to its owner GrunObject (example: 1L2W, A, V, or can even be simply a number like: 1, 8, 2.5 etc)
*		_relationQuantity:			the result value when the item's _relationship is parsed and evaluated (represents how much of the GrunObject the item deals with in spacial measurement units like m, m2 or m3 for examples)
*		_itemQuantityFormula:		factors/terms applied to the _relationQuantity. _relationQuantity is the LHS of the math expression, _itemQuantityFormula is the RHS. (Example: relationQuantity -> [50.0]m2 [/ 12.5]m2/mat of mesh <- quantityFormula to create the result expression: (50.0 / 12.5) which is the itemQuantity: 4.0 (mats of mesh))
*		_itemQuantity:				the result value when the relationQuantity & quantityFormula are parsed and evaluated.
*		_itemPrimaryLabourFormula:	factors/terms applied to the _itemQuantity. _itemQuantity is the LHS of the math expression, _itemPrimaryLabourFormula is the RHS. (Example: itemQuantity -> [4.0]mats of mesh [/ 1.5]mats laid/hour <- primaryLabourFormula to create the result expression: (4.0 / 1.5) which is the itemPrimaryLabour: 2.66 (hours))
*		_itemPrimaryLabour:			the result value when the itemQuantity & primaryLabourFormula are parsed and evaluated.
 */
class GrunItem
{
	public:
	static const GrunItem					INVALID;													// a static instance of an invlaid GrunItem
	std::function<void(GrunItem&)>			onValueChange;												// a function pointer that takes a reference to an instance of this class - this is used to recalculate GrunItem values onValueChange

	int										_libraryId					= -1;							// maps to the 'id' field in the UserInventory database table (this 'proves' the item's information has come from a valid source. if this value is -1 it means the object data is invalid)
	std::string 							_itemName					= "";							// required value on construction
	std::vector<RelationshipValues>			_itemCoreValues				= {};							// the core values in a GrunItem that must stay synced together

	// simple calculated members (these member values are calculated simply from the _itemQuantity value)
	double									_spatialTotalQuantity		= 0.0;							// the total 'Spatial Quantity' value for all relationships in the GrunItem's instance
	double									_itemTotalQuantity			= 0.0;							// the total 'Item Quantity' value for all relationships in the GrunItem's instance
	double									_itemTotalPrimaryLabour		= 0.0;
	double									_itemTotalQuantityRounded	= 0.0;
	double									_itemWasteFactor			= 0.0;
	double									_itemWasteAllowance			= 0.0;
	double									_itemItemizedProfitFactor	= 0.0;
	double									_itemItemizedProfit			= 0.0;
	
	// 'outsourced' members (they do not change based on the owning Grunobject's properties, usually stored in an inventory or database and passed in when the GrunItem is created)
	// overloaded ctr's will allow the dev to supply additional GrunItem values on instantination - this is for future development when GrunItem data will be sourced from a database or the network
	std::string 							_itemPrimaryLabourFormula	= "";
	std::string 							_itemQuantityFormula		= "";
	std::string								_itemCategory				= "";
	double									_itemRoundUpFactor			= 1.0;
	std::string								_itemSupplier				= "";
	std::string 							_itemSupplierSKU			= "";
	std::string								_itemSupplierDescription	= "";
	long long								_itemCostPerUnitCents		= 0;										// cost per GrunItem unit stored as signed 64 bit integer (long long)
	std::string 							_itemQuantityUnits			= "unit(s)";
	std::string								_itemPrimaryLabourUnits		= "hour(s)";
	bool									_hideFromClientView			= false;
	std::string 							_clientViewMessage			= "Install " + _itemName + " as required";
	std::chrono::system_clock::time_point	_itemLKGWUpdated{};
	std::chrono::system_clock::time_point	_itemLKGWCalculated{};

	// values that require outsourced and calculated values are defined below here
	

	/* Additional ideas for more GrunItem attributes that are beyond what exists in the Spreadsheet are: 
		- relationships to outside sources (databases) for relevant information such as:
			+ Material Safety Data (MSDS - if there's any safety data for this Item it can be included in the Construction Project's scope of documentation or made easily accessible digitally to all onsite workers)
			+ Safe Works Method Statement (SWMS - any relevant data about this item and its hazard ratings for its Primary Labour needs to be included into the Construction Project's specific SWMS. The SWMS can be made digitally available to all onsite workers)
			+ General, or In-House Specific Training information (The Grun user (the project manager/business owner) can link the item to an educational/informational URL that onsite workers can digitally access for training purposes)
			+ Schedule of In-House Owned Assets/Required Consumables for the GrunItem's installation (example: If you have a 400mm Shuttered Formwork GrunItem, the software could be told that a 400mm shutter is comprised of 400mm ply strips (lineal m), 2 * 4x2 LVLs (lineal m), 2 * 900 picket (per lineal m), 3 * tek screws (per lineal m), and 3 * 2 inch nails (per lineal m))
	*/

	// default ctr - assigns values to required fields
	GrunItem(	std::string name,
				int			id						= -1,
				std::string relationship			= "",
				std::string	relComment				= "",
				std::string quantityFormula			= "", 
				std::string units					= "unit(s)",
				std::string primaryLabourFormula	= ""
			);

	std::string	getCalculatedTimeString(const std::chrono::system_clock::time_point& member, const std::string& format = "%Y%m%d %H:%M:%S");
	int 		getLibraryId() const { return _libraryId; }
	int			getNumberOfRelationships() { return _itemCoreValues.size(); }
	
	/**
	 * @brief Add new RelationshipValues to GrunItem. Allows default values.
	 * @param	relationship			(std::string)	a relationship string for the owning GrunItem
	 * @param	relComment				(std::string)	a comment string specifically to describe the relationship (optional)
	 * @param	itemQuantity			(double)		the calculated item quantity based on the relationship
	 * @param	spatialQuantity			(double)		the calculated spatial quantity based on the relationship
	 * @param	spatialUnit				(SpatialExponentValue)	the calculated spatial unit based on the relationship
	 */
	size_t 		addRelationshipValues(	const std::string& 			relationship	= "", 
										const std::string& 			relComment		= "",
										const double& 				itemQuantity 	= 0.0, 
										const double& 				spatialQuantity	= 0.0, 
										const SpatialExponentValue& spatialUnit 	= SpatialExponentValue::None
								);

	size_t		rmRelationshipValues(	const size_t				index);
	
	GrunItem&	updateRelationshipValue(const size_t				index, 
										const std::string			relationship,
										const std::string			relComment);

	// Overload the equality operator
    bool operator==(const GrunItem& other) const 
	{	
        // We check the Library ID as the unique identifier
        return this->_libraryId == other._libraryId;
    }

	bool operator!=(const GrunItem& other) const { return !(*this == other); }
};
