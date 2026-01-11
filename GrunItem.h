#include <chrono>
#include <ctime>

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
/**
 * @brief CoreValues struct are the members in Grunitem that are all linked together via the Item's relationship.
 * @note Because we need to support Items having multiple relationships, this struct will help us keep all members that are influenced by the item's relationship in a neat unit.
 */
	struct CoreValues
	{
		std::string				relationship;
		std::string				relComment;
		bool					isCompoundRelationship;
		double					itemQuantity;
		double					spatialQuantity;
		SpatialExponentValue	spatialUnit;
	};

	public:
	std::string 							_itemName					= "";							// required value on construction
	std::vector<CoreValues>					_itemCoreValues				= {};							// the core values in a GrunItem that must stay synced together
	std::string								_relationship;												// the relationship will ultimately NOT be required on object creation, only the itemName is
	std::string								_comment					= "";							// a comment hte end user can put in for the item

	// interpretted member values (members that are derived when the item's _relationship is interpretted)
	// the members are listed in a rough locigal order of when they are calculated in the code
	bool									_isCompoundRelationship		= false;						// _isCompoundRelationship is determined by checking if the GrunItem's _calculatedSpatialUnit is smaller than its _itemQuantitySpatialUnit
	SpatialExponentValue					_spatialUnit				= SpatialExponentValue::None;	// the 'Spatial Unit' value (after interpretting and considering the entire Base Expression)
	double									_spatialQuantity			= 0.0;							// the 'Spatial Quantity' value
	SpatialExponentValue					_itemQuantitySpatialUnit	= SpatialExponentValue::None;	// _itemQuantitySpatialUnit is the SpatialExponentValue (None,Linear,Area,Volume) that is assigned to the GrunItem based on the GrunItem's _itemQuantityUnits value *IF* the _itemQuantityUnits are already of a spatial unit type (dev-note: this mostly only works if the _itemQuantityUnits are 'm', 'm2', 'm3' and these values are hard coded in GrunObject::mapUnitToSpatialExponent() which will need to be more flexible for locales in the future)
	double									_itemQuantity				= 0.0;

	// development members (these are mostly here for during development to check things)
	std::string								_baseExpression				= "";							// _baseExpression is the portion of the GrunItem's relationship string that is interpretted to result in the GrunItem's _calculatedSpatialUnit
	std::string								_baseExpressionIntprForSU	= "";							// the base expression interpretted for calculating the Spatial Unit
	std::string								_baseExpressionIntprNumeric	= "";							// the base expression after interpretation with numeric values in place of GrunObject Tokens and + in place of * operators
	std::string								_interprettedRelationship	= "";							// the interpretted relationship of the GrunItem for Item Qty calculation purposes
	SpatialExponentValue					_spatialAnchor				= SpatialExponentValue::None;	// the 'Spatial Anchor' value
	std::string								_spatialQuantityFormula		= "";							// the formula that is derived by converting the Base Expression's SHN into numeric math formula

	// simple calculated members (these member values are calculated simply from the _itemQuantity value)
	double									_itemPrimaryLabour			= 0.0;
	double									_itemRoundUpFactor			= 1.0;
	double									_itemQuantityRounded		= 0.0;
	double									_itemWasteFactor			= 0.0;
	double									_itemWasteAllowance			= 0.0;
	double									_itemItemizedProfitFactor	= 0.0;
	double									_itemItemizedProfit			= 0.0;
	
	// generally 'static' members (they do not change based on the owning Grunobject's prooperties, usually stored in an inventory or database and passed in when the GrunItem is created)
	// overloaded ctr's will allow the dev to supply additional GrunItem values on instantination - this is for future development when GrunItem data will be sourced from a database or the network
	std::string 							_itemPrimaryLabourFormula	= "";
	std::string 							_itemQuantityFormula		= "";
	std::string								_itemCategory				= "";
	std::string								_itemSupplier				= "";
	std::string 							_itemSupplierSKU			= "";
	std::string								_itemSupplierDescription	= "";
	double									_itemCostPerUnit			= 0.0;
	std::string 							_itemQuantityUnits			= "unit(s)";
	std::string								_itemPrimaryLabourUnits		= "hour(s)";
	bool									_hideFromClientView			= false;
	std::string 							_clientViewMessage			= "";
	std::chrono::system_clock::time_point	_itemLKGWUpdated{};
	std::chrono::system_clock::time_point	_itemLKGWCalculated{};
	// add more Item attributes as you need them through development

	/* Additional ideas for more GrunItem attributes that are beyond what exists in the Spreadsheet are: 
		- relationships to outside sources (databases) for relevant information such as:
			+ Material Safety Data (MSDS - if there's any safety data for this Item it can be included in the Construction Project's scope of documentation or made easily accessible digitally to all onsite workers)
			+ Safe Works Method Statement (SWMS - any relevant data about this item and its hazard ratings for its Primary Labour needs to be included into the Construction Project's specific SWMS. The SWMS can be made digitally available to all onsite workers)
			+ General, or In-House Specific Training information (The Grun user (the project manager/business owner) can link the item to an educational/informational URL that onsite workers can digitally access for training purposes)
			+ Schedule of In-House Owned Assets/Required Consumables for the GrunItem's installation (example: If you have a 400mm Shuttered Formwork GrunItem, the software could be told that a 400mm shutter is comprised of 400mm ply strips (lineal m), 2 * 4x2 LVLs (lineal m), 2 * 900 picket (per lineal m), 3 * tek screws (per lineal m), and 3 * 2 inch nails (per lineal m))
	*/

	// default ctr - assigns values to required fields
	GrunItem(	std::string name,
				std::string relationship = "",
				std::string quantityFormula = "", 
				std::string units = "unit(s)",
				std::string primaryLabourFormula = ""
			);

	std::string	getCalculatedTimeString(const std::chrono::system_clock::time_point& member, const std::string& format = "%Y%m%d %H:%M:%S");
	int			getNumberOfRelationships() { return _itemCoreValues.size(); }
	
	/**
	 * @brief Set GrunItem's core values. Allows default values.
	 * @param	relationship			(std::string)	a relationship string for the owning GrunItem
	 * @param	relComment				(std::string)	a comment string specifically to describe the relationship (optional)
	 * @param	isCompoundRelationship	(bool)			if the relationship is compound or not	
	 * @param	itemQuantity			(double)		the calculated item quantity based on the relationship
	 * @param	spatialQuantity			(double)		the calculated spatial quantity based on the relationship
	 * @param	spatialUnit				(SpatialExponentValue)	the calculated spatial unit based on the relationship
	 */
	void	setCoreValues(	std::string relationship			= "", 
							std::string relComment				= "",
							bool isCompoundRelationship			= false, 
							double itemQuantity 				= 0.0, 
							double spatialQuantity				= 0.0, 
							SpatialExponentValue spatialUnit 	= SpatialExponentValue::None
						);
};
