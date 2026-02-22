#include "GrunStage.h"

/**
 * @brief Construct a new Grun Stage:: Grun Stage object
 * 
 * @param stageName	(std::string)	- the GrunStage instance's Stage name
 */
GrunStage::GrunStage(const std::string& stageName) : m_name(stageName)
{
	// add any involved construction code here as needed
}

/**
 * @brief adds a GrunObject to the GrunStage instance
 * 
 * @param gObject (GrunObject)	- the GrunObject to remove
 * @return size_t Returns the number of GrunObjects now in the m_objects member
 */
size_t GrunStage::addGrunObject(const GrunObject& gObject)
{
	m_objects.push_back(gObject);
	return m_objects.size();
}

StageResults GrunStage::calculateMVPAreaRate(const bool rounded)
{
	StageResults results;
	for (GrunObject& obj : m_objects)
	{
		results.totalCostCents += obj.getTotalCostOfObject(rounded);
	}
	return StageResults();
}

/**
 * @brief Calculates and updates the m_rateCost member's value (the total cost to construct this Stage)
 * 
 * @return (size_t)			- the rateCost as an unsigned integer (no decimal places)
 */
size_t GrunStage::calculateRateCost()
{
	size_t total = 0;
	for (const auto& object : m_objects)
	{
		object;
	}
	return size_t();
}

size_t GrunStage::createGrunObject(const std::string &shapeType, const std::string &name, double x, double y, double z, const std::string &areaType, const std::string &stage)
{
	// use GrunObject ctr to create the new GrunObject - we can rely on default values if no values are given in the arguments
	GrunObject	newGrunObject(shapeType, name, x, y, z, areaType, stage);
	return this->addGrunObject(newGrunObject);
}

std::optional<std::reference_wrapper<GrunObject>> GrunStage::getGrunObject(const size_t index)
{
	// zero & out of bounds check
	if (index >= m_objects.size()) { return std::nullopt; }
	return m_objects[index];
}

const std::string GrunStage::getName()
{
		return m_name;
}

/**
 * @brief Remove a GrunObject from the GrunStage
 *
 * @param objIndex (int)	- the inde of the GrunObject to remove
 * @return size_t Returns the remaining number of GrunObjects in m_objects
 */
size_t GrunStage::rmGrunObject(const int& objIndex)
{
	m_objects.erase(m_objects.begin() + objIndex);
	return m_objects.size();
}

/**
 * @brief Set the GrunObject::m_name member
 * 
 * @param name (std::string)	- The new name
 */
void GrunStage::setName(const std::string& name)
{
	m_name = name;
	return;
}

/**
 * @brief Reutrns how many GrunObjects exist in this GrunStage instance
 * 
 * @return size_t	- the number of GrunItems that are housed in the instance's m_objects vector
 */
size_t GrunStage::size()
{
	return m_objects.size();
}

/**
 * @brief set the Stage's m_rateUnit value using an integer (only 0-3 are valid values)
 * 
 * @param spatialUnitAsInt	(size_t)	- integer value that represents the desired SpatialExponentValue for the Stage's rateUnit
 * @return (bool)						- returns true if the input argument was valid (0-3) otherwise returns false to indicate we set to 'None' due to a data input problem
 */
bool GrunStage::setRateUnit(const size_t &spatialUnitAsInt)
{
	if (spatialUnitAsInt <= 3)
	{
		m_rateUnit = static_cast<SpatialExponentValue>(spatialUnitAsInt);
		return true;
	}
	else
	{
		m_rateUnit = SpatialExponentValue::None;
		return false;
	}
}