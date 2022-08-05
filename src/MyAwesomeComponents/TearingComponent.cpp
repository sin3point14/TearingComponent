#include "TearingComponent.inl"
#include <sofa/defaulttype/VecTypes.h>

namespace sofa::component::controller
{
	int TearingControllerClass = core::RegisterObject("Adds tearing. kek")
		.add< TearingComponent<sofa::defaulttype::Vec3Types> >();

template class MYAWESOMECOMPONENTS_API TearingComponent<sofa::defaulttype::Vec3Types>;

}
