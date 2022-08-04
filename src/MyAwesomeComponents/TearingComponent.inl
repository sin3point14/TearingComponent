#include "TearingComponent.h"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ForceField.h>

namespace sofa::component::controller
{
	template <class DataTypes>
	TearingComponent<DataTypes>::TearingComponent()
		: d_yieldStress(initData(&d_yieldStress, (Real)10.0f, "yieldStress", "Yield Stress"))
		, m_triangleFEMForceField(nullptr)
	{
		this->f_listening.setValue(true);
	}

	template <class DataTypes>
	TearingComponent<DataTypes>::~TearingComponent()
	{
	}


	template <class DataTypes>
	void TearingComponent<DataTypes>::init()
	{

		// We look for a CollisionModel identified with the CarvingSurface Tag.

		m_triangleFEMForceField = getContext()->get<sofa::component::forcefield::TriangleFEMForceField<DataTypes>>();


		if (m_triangleFEMForceField == nullptr)
		{
			msg_error() << "No TriangularFEMForceField component found";
			sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
			return;
		}
	}


	template <class DataTypes>
	void TearingComponent<DataTypes>::reset()
	{
	}


	template <class DataTypes>
	void TearingComponent<DataTypes>::doTear()
	{
	}

}
