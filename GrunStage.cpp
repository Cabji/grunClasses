#include "GrunStage.h"

/**
 * @brief Construct a new Grun Stage:: Grun Stage object
 * 
 * @param stageName	(std::string)	- the GrunStage instance's Stage name
 */
GrunStage::GrunStage(const std::string& stageName) : m_name(stageName)
{

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

/**
 * @brief Remove a GrunObject from the GrunStage
 * 
 * @param objIndex (int)	- the inde of the GrunObject to remove
 * @return size_t Returns the remaining number of GrunObjects in m_objects
 */
size_t GrunStage::rmGrunObject(const int& objIndex)
{
	m_objects.erase(objIndex);
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
