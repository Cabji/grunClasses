# Importing CSV data to MVP UserInventory Table

## Ensure your import CSV data matches the fiels in the UserInventory table

When you want to import CSV data into the MVP UserInventory database table, you need to make the import data's format match the field types in the UserInventory table schema.

### The UserInventory table schema

You should check directly in the source code to ensure you have the latest design of the table's schema as it may change over time. The table's schema can be found in the Grun project source code: DatabaseManager::initializeSchema() in DatabaseManager.cpp.

This file was written on 20260208 and the schema was, at time of writing, the following: 

				  "CREATE TABLE IF NOT EXISTS UserInventory ("
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

## Format of a record in the import CSV file

The format for a record in the import data is: 
id,itemName,itemCategory,units,itemQtyFormula,itemCost,roundUpFactor,primaryLabourFormula,primaryLabourUnits,hideFromClientView,itemClientViewMessage,lkgwUpdateTime

A sample row in the import CSV is:
NULL,"End User Seen Item Name","End User Seen - Category","unit(s)","",4850,1.0,"","hour(s)",0,"Install Item Name as required in proposed project",NULL

where NULL values will be handled by sqlite using default values

## Importing the CSV Data

When you import CSV data to the UserInventory table, you should use the command below, copy and pasted into a terminal. I am using the msys2 UCRT terminal to do this on Windows, so this should also work under Linux environment.

```
sqlite3 YourSQLite3DatabaseFilename.db <<EOF
.mode csv
.nullvalue NULL
CREATE TEMP TABLE staging(id, name, cat, unit, form, cost, rnd, lab_form, lab_unit, hide, msg, ts);
.import your_import_csv_filename.csv staging

-- The WHERE clause filters out that pesky blank 11th line
INSERT INTO UserInventory (item_name, item_category, item_qty_unit, item_qty_formula, item_cost_per_unit_cents, item_round_up_factor, item_primary_labour_formula, item_primary_labour_units, default_hide_from_client_view, item_client_view_message)
SELECT name, cat, unit, form, cost, rnd, lab_form, lab_unit, hide, msg 
FROM staging 
WHERE name IS NOT NULL AND name != '';

DROP TABLE staging;
SELECT id, item_name FROM UserInventory;
EOF
```
This 'staging' table method was required because SQLite would always complain about the NULL values in the import data, even though those fields are meant to have default values if they are not provided. Gemmy advised that a staging table is the "industry standard" way to do this.