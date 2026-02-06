#include "DatabaseManager.h"

DatabaseManager::DatabaseManager(const std::string& dbPath)
	: m_db(dbPath, SQLite::OPEN_READWRITE | SQLite::OPEN_CREATE)
{
}

void DatabaseManager::initializeSchema() 
{
    SQLite::Transaction transaction(m_db);

    // 1. The Core Product Catalog
    m_db.exec("CREATE TABLE IF NOT EXISTS MasterProducts ("
              "id INTEGER PRIMARY KEY AUTOINCREMENT,"
              "name TEXT UNIQUE NOT NULL,"
              "category TEXT,"
              "base_unit TEXT NOT NULL"
              ");");

    // 2. The Supplier Price Matrix (The "Marketplace")
    m_db.exec("CREATE TABLE IF NOT EXISTS SupplierOffers ("
              "id INTEGER PRIMARY KEY AUTOINCREMENT,"
              "product_id INTEGER NOT NULL,"
              "supplier_name TEXT NOT NULL,"
              "price_level TEXT DEFAULT 'Retail'," // Retail, Trade, Account
              "cost_cents INTEGER NOT NULL,"
              "FOREIGN KEY(product_id) REFERENCES MasterProducts(id),"
              "UNIQUE(product_id, supplier_name, price_level)"
              ");");

    // 3. User Customizations (The "Aliases")
    m_db.exec("CREATE TABLE IF NOT EXISTS UserInventory ("
              "id INTEGER PRIMARY KEY AUTOINCREMENT,"
              "master_product_id INTEGER NOT NULL,"
              "alias_name TEXT UNIQUE,"
              "preferred_offer_id INTEGER,"
              "manual_price_override_cents INTEGER," // NULL if using Supplier price
              "FOREIGN KEY(master_product_id) REFERENCES MasterProducts(id),"
              "FOREIGN KEY(preferred_offer_id) REFERENCES SupplierOffers(id)"
              ");");

    transaction.commit();
}