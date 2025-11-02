#include <chrono>
#include <ctime>

#ifndef CLASS_NAME
#define CLASS_NAME "GrunItem"
#endif

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
	std::string 							_itemName					= "";		// required value on construction
	std::string 							_relationship				= "";		// required value on construction
	double									_relationQuantity			= 0.0;
	std::string 							_itemQuantityFormula		= "";
	double									_itemQuantity				= 0.0;
	std::string 							_itemPrimaryLabourFormula	= "";
	double									_itemPrimaryLabour			= 0.0;
	// all values above, excluding _itemQuantity and _itemPrimaryLabour are required on GrunItem instantination as a minimum for the default construction

	// overloaded ctr's will allow the dev to supply additional GrunItem values on instantination - this is for future development when GrunItem data will be sourced from a database or the network
	double									_itemRoundUpFactor			= 1.0;
	double									_itemQuantityRounded		= 0.0;
	double									_itemWasteFactor			= 0.0;
	double									_itemWasteAllowance			= 0.0;
	double									_itemItemizedProfitFactor	= 0.0;
	double									_itemItemizedProfit			= 0.0;
	
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
	GrunItem(std::string name, std::string relationship, std::string quantityFormula = "", std::string units = "unit(s)", std::string primaryLabourFormula = "")
		: _itemName(name), _relationship(relationship), _itemQuantityFormula(quantityFormula), _itemQuantityUnits(units), _itemPrimaryLabourFormula(primaryLabourFormula)
	{
		// add any needed calculations in here, example: use _relationship to determine value that allows quantifying the GrunItem's _itemQuantity
	}

	/**
	 * @brief Converts a GrunItem's time-typed member to user-friendly date/time string, returned in optional format
	 * @param member 	- the time-typed member in the GrunItem we want to retrieve (required)
	 * @param format	- std::string that defines the format the time should be shown in (default: "%Y%m%d %H:%M:%S")
	 * @return std::string	- The formatted datetime string, NULL" if the item's attributes have never been calculated, or error msg if an error is encountered
	 */
	std::string	getCalculatedTimeString(const std::chrono::system_clock::time_point& member, const std::string& format = "%Y%m%d %H:%M:%S") 
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
};
