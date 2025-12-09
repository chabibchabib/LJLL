static TypeOfFE_PkLagrange PKLagrange(k); // Obj of the class TypeOfFE_PkLagrange 
namespace Fem2D {
    ...
    static void init() {

        AddNewFE("PKLagrange", &PKLagrange); 
        static ListOfTFE FE_PK("PKLagrange", &PKLagrange); //add Pk to the list of Common FE
    }
} // namespace Fem2D
LOADFUNC(Fem2D::init);