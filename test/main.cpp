#include <cstdlib>
#include <filesystem>
#include <format>
#include <iostream>
#include <print>
#include <ranges>
#include <string>
#include <vector>
#include "GrunProject.h"
#include "DatabaseManager.h"
#include "InventoryManager.h"


int main ()
{
	try 
	{
		std::string dbFilename = get_app_data_path().string().append("/GrunInventory.db");	// get the full path to the database file
        SQLite::Database db(dbFilename, SQLite::OPEN_READONLY);								// create a database connection instance
        InventoryManager inventory(db);														// create an InventoryManager instance

		GrunProject	projectInstance("Test Project");										// create a GrunProject instance
		GrunStage& stageFootings = projectInstance.createStage("Footings");					// create a GrunStage in the GrunProject instance
		stageFootings.createGrunObject("rectangle","SF1",50.0,0.4,0.4);						// create a GrunObject in the Footings stage object (instance)
		stageFootings.createGrunObject("rectangle","PF1",1.0,1.0,0.3);						// create a GrunObject for a Pad Footing in the stage object
		auto wrappedObj = stageFootings.getGrunObject(0);
		if (wrappedObj)
		{
			// now, the GrunProject instance has some sort of access to an hourly rate
			size_t hourlyRate = projectInstance.projectLabourRates.tiers[0].hourlyRate;
			GrunObject& sf1Obj = wrappedObj.value();
			// user InventoryManager object to fuzzy search and find best matches for your items in the database, then add the found GrunItem using GrunObject::addItem(GrunItem&)
			sf1Obj.addGrunItem("Excavator",1000,"8","","","hour(s)","/2");
			sf1Obj.addGrunItem(inventory.getBestMatch("T M 11 3")).updateRelationshipValue(0,"1L","SF1 Steel");
			sf1Obj.addGrunItem(inventory.getBestMatch("chair multi")).updateRelationshipValue(0,"1L@0.8","SF1 Trench Chairs");
			sf1Obj.addGrunItem(inventory.getBestMatch("CONC 20/20 footings")).updateRelationshipValue(0,"V","SF1 Concrete");
			GrunItem& starters = sf1Obj.addGrunItem(inventory.getBestMatch("starter bar 1200"));
			starters.updateRelationshipValue(0,"1L@0.8","SF1 Steel - Slab Ties @0.8");
			starters.addRelationshipValues("0.75L@0.4","SF1 Steel - BW1 Starters @0.4");
			starters.addRelationshipValues("2 * 4","SF1 Steel - Pier Extra Bars (2 piers, 4 bars per pier)");
			// when relationships are added to GrunItems, the calculations are done at the GrunItem scope to determine how much of each item is required
			// we need to Calculate the totals at the GrunObject level
			//std::println("{}",sf1Obj.getGrunItemListInfoAsString("%Y%m%d",32));
		}
		else
		{
			std::println("Failed to find an object in m_objects[0]. Did you acutally add any objects to the stage?");
		}
		wrappedObj = stageFootings.getGrunObject(1);
		if (wrappedObj)
		{
			GrunObject& pf1Obj = wrappedObj.value();
			// user InventoryManager object to fuzzy search and find best matches for your items in the database, then add the found GrunItem using GrunObject::addItem(GrunItem&)
			pf1Obj.addGrunItem(inventory.getBestMatch("NF82")).updateRelationshipValue(0,"5*2A","PF1 Mesh Top & Bottom");
			pf1Obj.addGrunItem(inventory.getBestMatch("5065")).updateRelationshipValue(0,"5*A","PF1 Bottom Mesh Chairs 50/65 SOG");
			pf1Obj.addGrunItem(inventory.getBestMatch("CONC 20/20 footings")).updateRelationshipValue(0,"5*V","PF1 Concrete");
		}
		else
		{
			std::println("Failed to find an object in m_objects[1]. Did you acutally add any objects to the stage?");
		}
		std::println("Total Cost for Stage '{}': $ {}", stageFootings.getName(),stageFootings.calculateRateCost() / 100);

    }
    catch (const std::exception& e) {
        std::println(std::cerr, "CRITICAL ERROR: {}", e.what());
        return 1;
    }

    return 0;

/* test the Inventory Manager
	try {
        // 1. Initialize Database
		std::string dbFilename = get_app_data_path().string().append("/GrunInventory.db");
        SQLite::Database db(dbFilename, SQLite::OPEN_READONLY);
        InventoryManager inventory(db);

        std::string searchTerm;
        std::println("--- Grun Bare Basics Inventory Test (C++23) ---");

        while (true) {
            std::print("\nEnter item to search (or 'exit' to quit): ");
            if (!std::getline(std::cin, searchTerm) || searchTerm == "exit") break;

            // 2. Perform the Bare Basics Search
            auto results = inventory.search(searchTerm);

            if (results.empty()) {
                std::println("No items found matching '{}'", searchTerm);
                continue;
            }

            // 3. Display Results and verify Hydration
            std::println("Found {} items:", results.size());
            
            for (const auto& item : results) {
                std::println("{:-<42}", ""); // Horizontal separator
                std::println("Name:     {}", item._itemName);
                std::println("Category: {}", item._itemCategory);
                
                // Formatting currency to 2 decimal places using format specifiers
                std::println("Cost:     ${:.2f} per {}", 
                             (item._itemCostPerUnitCents / 100.0), 
                             item._itemQuantityUnits);
                
                std::println("Formula:  {}", item._itemQuantityFormula);
                std::println("Labour:   {} {}", 
                             item._itemPrimaryLabourFormula, 
                             item._itemPrimaryLabourUnits);
            }
        }
    }
    catch (const std::exception& e) {
        std::println(std::cerr, "CRITICAL ERROR: {}", e.what());
        return 1;
    }

    return 0;
*/

/*
// DatabaseManager testing
	fs::path	configFolder	= get_app_data_path();
	if (configFolder.empty()) 
	{
		std::println("Fatal Error: Could not create program config folder in user's home folder, so quitting.");
		return 1;
	}

	std::string	dbPath	= (configFolder / "GrunInventory.db").string();
	// create a database instance
	DatabaseManager db(dbPath);
	db.initializeSchema();

// Grun Classes testing	
	// create a vector to store the stages (we will make a Project class later on)
	std::vector<GrunStage> project;
	// create a GrunStage for the Footings stage of this project
	GrunStage	stageFootings("Footings");

	// imagine we want to build a 9.0x6.0 garage slab, with a strip footing below it

	// create a GrunObject (3D Solid) in the Footing stage object that represents the strip footing in 3D space
	stageFootings.createGrunObject("rectangle", "SF1", 30.0, 0.3, 0.4, "horizontal");
	// temporarily get a reference the 3D Solid
	auto objWrap = stageFootings.getGrunObject(0);
	if (objWrap) 
	{ 
		// add GrunItems to the 3D Solid via the reference we made
		auto& object = objWrap->get(); 
		object.addGrunItem("Excavator","8","","","hour(s)","*0.5");
		object.addGrunItem("Trench Mesh - 3 Bar TM11","1L","","/5.4","length(s)","/6");
		object.addGrunItem("Chairs - Trench Mesh","1L@0.8","","/25","bag(s) of 25","/4");
		object.addGrunItem("Slab Ties - N12","1L@0.4","1.2x0.25 L Bar N12","","bar(s)","/30");
		object.addGrunItem("Delivery - Steel","1","To site address","","delivery(ies)","*0.5");
		object.addGrunItem("Concrete - 20/20 (Footings)","V","","","m3","*1.55");
		// now the program knows what's needed to construct the SF1 Footing
	}

	// create a GrunStage for the Slab that will be built on top of the footing
	GrunStage stageSlab("Slab");
	// create a GrunObject (3D Solid) in the Slab stage object that represents the slab in 3D space
	stageSlab.createGrunObject("rectangle", "Slab", 9.0, 6.0, 0.1, "horizontal");
	// temporarily get a reference the 3D Solid
	objWrap = stageSlab.getGrunObject(0);
	if (objWrap) 
	{ 
		// add GrunItems to the 3D Solid via the reference we made
		auto& object = objWrap->get(); 
		object.addGrunItem("Delivery - Subgrade","1","To site address","","delivery(ies)","*0.5");
		object.addGrunItem("Subgrade - Fines","0.05A","","","m3","*1");
		object.addGrunItem("Formwork - 300 Shuttered","2L2W","","","m","/2.46");
		object.addGrunItem("Viscrine 4mx50m","A","","/200","roll(s)","*1");
		// dev-note: new idea - reference to another item in the project in the relationship to get the quantity of that item
		object.addGrunItem("Secondary Labour","&SF1::Slab Ties N12","Bend slab ties into slab","","","/30");
		object.addGrunItem("Mesh SL82 (6.0x2.4)","A","","/12.5","mat(s)","/1.5");
		object.addGrunItem("Chairs SOG 5065","A","","/(0.6*0.6)/100","bag(s) of 100","*0.5");
		object.addGrunItem("Delivery - Steel","1","To site address","","delivery(ies)","*0.5");
		object.addGrunItem("Concrete - 20/20","V","","","m3","*3.2");
		object.addGrunItem("Strip Formwork","2L2W","","","m","/6");
		// now the program knows what's needed to construct the Slab
	}

	// push the GrunStages into our project vector so we can loop through them for output
	project.push_back(stageFootings);
	project.push_back(stageSlab);

	std::println("Project has {} stages",project.size());
	for (auto& stage : project)
	{
		std::println("Stage: {}",stage.getName());
		for (auto i = 0; i < stage.size(); ++i)
		{
			// get the current GrunObject
			auto wrappedObject = stage.getGrunObject(i);
			if (!wrappedObject) { continue; }
			auto& currentObject = wrappedObject->get();
			std::println("GrunObject [{}] Item List information:\n{}", currentObject.getObjectName(), currentObject.getGrunItemListInfoAsString("%Y%m%d %H%M%S"));
		}
	}
*/

/* this comment block demonstrates how to use the search methods to find relationships and their owning GrunItem locations in GrunObject::m_items
	// search all items, for a relationship
	std::string comSearch									= "";
	std::string relSearch									= "?L?W";
	// create a vector with all the indices from the slab's (GrunObject's) m_items property, using modern C++23 views and ranges
	std::vector<size_t>	searchIndices						= std::views::iota(0uz, slab.getTotalOfGrunItems()) | std::ranges::to<std::vector>();
	std::vector<RelationshipSearchResult>	searchResults	= slab.findRelationshipByStrings(searchIndices,relSearch,comSearch,false);

	std::print("Search results for Rel:'{}' and Com:'{}':",relSearch,comSearch);
	if (!searchResults.empty()) 
	{ 
		for (const auto result : searchResults)
			std::print(" {}", result);
	}

	// std::print("Searched m_items{}. Searched for Relationship '{}' and Comment '{}'. Results were: {}", searchIndices, relSearch, comSearch, searchResults | std::views::all);
	if (searchResults.empty()) { std::println(" no relationship/comment matches found."); }
	std::println();
*/
	
	// std::println("GrunObject's details: Name: {}\n\tLength: {}\tWidth: {}\tDepth: {}\tArea: {}\tVolume: {}",slab.getObjectName(),slab.getObjectProperty("length"),slab.getObjectProperty("width"),slab.getObjectProperty("depth"),slab.getObjectProperty("area"),slab.getObjectProperty("volume"));
	// std::println("GrunObject [{}] Item List information:\n{}", slab.getObjectName(), slab.getGrunItemListInfoAsString("%Y%m%d %H%M%S"));
	// slab.calculateGrunObjectTotals();
	// std::println("GrunObject's Totals Data {}", slab.getGrunObjectTotalsInfoAsString());

	// std::println("Object's Item Total Data, as fetched by GrunObject::getItemQtyTotals()");
	// std::vector<ItemAndTotal> objTotals = slab.getItemQtyTotals(true);
	// for (const auto& itemTotal : objTotals)
	// {
	// 	std::println("{}",itemTotal);
	// }
}