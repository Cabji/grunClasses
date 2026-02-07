#include "DatabaseManager.h"


DatabaseManager::DatabaseManager(const std::string& dbPath)
    : m_db(dbPath, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE)
{
}

void DatabaseManager::initializeSchema() 
{
    // 0. Ensure we have a version tracking table
    m_db.exec("CREATE TABLE IF NOT EXISTS SchemaVersion (version INTEGER PRIMARY KEY);");
    
    // Get current version (default to 0 if table is empty)
    int currentVersion = 0;
    SQLite::Statement query(m_db, "SELECT version FROM SchemaVersion");
    if (query.executeStep()) {
        currentVersion = query.getColumn(0).getInt();
    } else {
        m_db.exec("INSERT INTO SchemaVersion (version) VALUES (0)");
    }

    SQLite::Transaction transaction(m_db);

    // --- Version 1: Initial Tables ---
    if (currentVersion < 1) {
        // The ONLY table needed for the Standalone MVP
    	m_db.exec("CREATE TABLE IF NOT EXISTS UserInventory ("
        	      "id INTEGER PRIMARY KEY AUTOINCREMENT,"							// primary key for database use
            	  "item_name TEXT UNIQUE NOT NULL,"									// user-friendly name of the item
				  "item_category TEXT,"												// category for the item
				  "item_qty_unit TEXT NOT NULL,"									// m, m2, m3, each
				  "item_qty_formula TEXT DEFAULT '',"								// the RHS of formula to convert from SpatialQty to ItemQty (e.g., '/ 12.5')
				  "item_cost_per_unit_cents INTEGER DEFAULT 0,"						// item cost per unit in cents (no decimals!)
				  "item_round_up_factor REAL DEFAULT 1.0,"							// the round up factor (0.2 for concrete, 1 for bag of chairs etc.)
				  "item_primary_labour_formula TEXT DEFAULT '',"					// the RHS of formula to convert from ItemQty to Primary Labour in labour units
				  "item_primary_labour_units TEXT DEFAULT 'hour(s)',"				// the units for the Primary Labour value
				  "default_hide_from_client_view INTEGER DEFAULT 0,"				// do we hide the item from the Client View by default?
				  "item_client_view_message TEXT DEFAULT '',"						// the message to show on the client view for this Item in the Job Logistics area
				  "lkgw_item_info_updated INTEGER DEFAULT (strftime('%s','now'))"	// timestamp in unix format of when the item's database record was last updated
				  ");");

        currentVersion = 1;
    }

    // --- Version 2: if we do updates we can do schema additions/adjustments here
    // if (currentVersion < 2) {

    //     currentVersion = 2;
    // }

    // Finalize the update
    SQLite::Statement updateVer(m_db, "UPDATE SchemaVersion SET version = ?");
    updateVer.bind(1, currentVersion);
    updateVer.exec();

    transaction.commit();
}