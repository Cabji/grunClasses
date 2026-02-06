#pragma once
#include <SQLiteCpp/SQLiteCpp.h>
#include <string>
#include <vector>

class DatabaseManager
{
	private:
	SQLite::Database						m_db;

	public:
					DatabaseManager(const std::string& dbPath);			// class ctr - open the database file, creates it if it doesn't exist
	void 			initializeSchema();									// create the initial table structure in the database
	void			executeCommand(const std::string& query);			// utility method to execute simple sql query
};