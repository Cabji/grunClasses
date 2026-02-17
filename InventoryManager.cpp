#include "InventoryManager.h"
#include "cabjiHelpers.h"

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
	item._itemCostPerUnitCents			= query.getColumn("item_cost_per_unit_cents").getInt64();
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
std::vector<GrunItem> InventoryManager::search(const std::string& nameQuery, const bool useFuzzySearch)
{
	std::vector<GrunItem>	results;
	std::string				fuzzySearchString = nameQuery;
	std::replace(fuzzySearchString.begin(), fuzzySearchString.end(), ' ', '%');

	SQLite::Statement query(m_db, "SELECT * FROM UserInventory WHERE item_name LIKE ? ORDER BY item_name ASC");
	(useFuzzySearch) ? query.bind(1, "%" + fuzzySearchString + "%") : query.bind(1, "%" + nameQuery + "%");
	

	while (query.executeStep())
	{
		GrunItem item(query.getColumn("item_name").getText(),query.getColumn("id").getInt64());
		hydrate(item,query);
		results.push_back(item);
	}
	return results;
}

GrunItem InventoryManager::getBestMatch(const std::string &nameQuery, const bool useFuzzySearch)
{
	// get the GrunItems that match the string via a fuzzy search
	std::vector<GrunItem> candidates = this->search(nameQuery, useFuzzySearch);
	
	// zero-check
	if (candidates.empty()) {
        // Log error or return a default "Not Found" item
        return GrunItem("Item Not Found", 0);
    }

	// pointer to GrunItem that is the best match
    GrunItem* bestMatch = nullptr;
    double highestScore = -1.0;

	// loop the cnadidate GrunItems
    for (auto& item : candidates) 
	{
		// the current GrunItem's _itemName
        std::string targetName = item._itemName;

        // calculate scores using different algorithms
        double tokenScore = cabji::get_token_score(nameQuery, targetName);
        double containScore = cabji::get_containment_score(nameQuery, targetName);
        double typoScore = cabji::get_levenshtein_score(nameQuery, targetName);

        // --- THE HYBRID WEIGHTING ---
        // Token match is most important for steel (60%)
        // Containment is great for "N12" matching "N12 REBAR" (30%)
        // Levenshtein catches "trech" vs "trench" (10%)
        double finalScore = (tokenScore * 0.6) + (containScore * 0.3) + (typoScore * 0.1);

        if (finalScore > highestScore) 
		{
            highestScore = finalScore;
            bestMatch = &item;
        }
    }

    return *bestMatch;
}
