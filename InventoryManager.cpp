#include "InventoryManager.h"

/**
 * @brief Hydrates a GrunItem with data from a result from the UserInventory
 * 
 * @param item 			(GrunItem&)				- reference to the GrunObject to hydrate
 * @param query 		(SQLite::Statement&)	- the query to fetch the item from the UserInventory table
 */
void InventoryManager::hydrate(GrunItem& item, SQLite::Statement& query)
{
	item._itemCategory					= query.getColumn("item_category").getText();
	item._itemQuantityUnits				= query.getColumn("item_qty_unit").getText();
	item._itemQuantityFormula			= query.getColumn("item_qty_formula").getText();
	item._itemCostPerUnit				= query.getColumn("item_cost_per_unit_cents").getInt64();
	item._itemRoundUpFactor				= query.getColumn("item_round_up_factor").getDouble();
	item._itemPrimaryLabourFormula		= query.getColumn("item_primary_labour_formula").getText();
	item._itemPrimaryLabourUnits		= query.getColumn("item_primary_labour_units").getText();
	item._hideFromClientView			= query.getColumn("default_hide_from_client_view").getInt() > 0;
	item._clientViewMessage				= query.getColumn("item_client_view_message").getText();
}

/**
 * @brief Search for a GrunItem's data in the UserInventory table ONLY (hard coded) by item_name field
 * 
 * @param nameQuery 			(std::string)	- the item name to search for (this is a fuzzy search using SQL's LIKE surrounded by % chars)
 * @return std::vector<GrunItem> 				- a vector of GrunItems that were found and hydrated
 */
std::vector<GrunItem> InventoryManager::search(const std::string& nameQuery)
{
	std::vector<GrunItem>	results;

	SQLite::Statement query(m_db, "SELECT * FROM UserInventory WHERE item_name LIKE ? ORDER BY item_name ASC");
	query.bind(1, "%" + nameQuery + "%");

	while (query.executeStep())
	{
		GrunItem item(query.getColumn("item_name").getText());
		hydrate(item,query);
		results.push_back(item);
	}
	return results;
}