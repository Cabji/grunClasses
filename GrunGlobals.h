#pragma once
#include <string>
#include <vector>

namespace fs = std::filesystem;

const size_t	DEFAULT_HOURLY_RATE	= 7500;

/**
 * @brief Gets the app data path object which is used to find where the program will store its needed data files
 * 
 * @return fs::path 
 */
fs::path	get_app_data_path() 
{
	fs::path	home;
#ifdef _WIN32
	// windows based
	const char*	drive	= std::getenv("HOMEDRIVE");
	const char*	path	= std::getenv("HOMEPATH");
	if (drive && path)
		home = fs::path(drive) / path;
	else
	{
		const char*	userProfile	= std::getenv("USERPROFILE");
		if (userProfile)
			home = fs::path(userProfile);
	}
#else
	// unix based
	const char*	homeEnv	= std::getenv("HOME");
	if (homeEnv)
		home = fs;:path(homeEnv);
#endif

	// append the app's config folder
	fs::path	appPath	= home / ".grun";
	// create directory if it doesn't exist
	if (!home.empty() && !fs::exists(appPath))
		fs::create_directories(appPath);

		return appPath;
}


/**
 * @brief A simple Labour Tier data strcture
 * 
 */
struct LabourTier
{
	std::string		name;
	size_t			hourlyRate;
};

/**
 * @brief The LabourRates struct lets us define an artbitrary number of LabourTiers housed in a vector
 * 
 */
struct LabourRates
{
	std::vector<LabourTier>	tiers;

	// ctr
	LabourRates()
	{
		tiers.push_back({"Default Hourly Rate", DEFAULT_HOURLY_RATE});
	}

	// helper to get Default Labour Tier object
	const LabourTier& getDefaultTier()
	{
		return tiers[0];
	}
};

class GrunGlobals
{
public:

	/**
	 * @brief The core of the Singleton
	 * Access via: GrunGlobals::getInstance().labourRates
	 * @return GrunGlobals& 
	 */
	static GrunGlobals& getInstance()
	{
		static GrunGlobals instance;
		return instance;
	}

	LabourRates	labourRates;

	// add other global defaults here later

private:
	// private ctr to prevent instantination elsewhere
	GrunGlobals() {}

	// delete copy constructor & assignment operator to prevent duplicate
	GrunGlobals(const GrunGlobals&) = delete;
	void operator=(const GrunGlobals&) = delete;
};