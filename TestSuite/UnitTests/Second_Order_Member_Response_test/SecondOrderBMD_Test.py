
import pytest
import numpy as np

from config import config
from StructuralElements import Node, Member
from Loads import NeumanBC
from FirstOrderResponse import FirstOrderMemberResponse
from SecondOrderResponse import  SecondOrderMemberResponse
from FiniteElementDivisor import divide_into_finite_elements


@pytest.fixture
def setup_model():
    config.set_FEDivision(11)
    PointsT = [
            Node(Node_Number=1, xcoordinate=0, ycoordinate=0, Support_Condition="Hinged Support"),
            Node(Node_Number=2, xcoordinate=0, ycoordinate=5, Support_Condition="Rigid Joint"),
            Node(Node_Number=3, xcoordinate=5, ycoordinate=5, Support_Condition="Hinged Support")
            ]


    MembersT = [
            Member(Beam_Number=1, Start_Node=PointsT[0], End_Node=PointsT[1], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
            Member(Beam_Number=2, Start_Node=PointsT[1], End_Node=PointsT[2], Area=0.09, Youngs_Modulus=200000000, Moment_of_Inertia=0.000675),
            ] # square cross section - 0.3 x 0.3, units N, m


    LoadsT = [
            NeumanBC(type="PL", Magnitude=-100000, Distance1=2.5, AssignedTo="Member 2", Members = MembersT)
            ] 

    PointsT, MembersT, LoadsT = divide_into_finite_elements(PointsT, MembersT, LoadsT, 10)

    FirstOrderResponse = FirstOrderMemberResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)
    SecondOrderResponse = SecondOrderMemberResponse(Points = PointsT, Members = MembersT, Loads = LoadsT)
    
    return FirstOrderResponse, SecondOrderResponse


def test_BendingMomentAll(setup_model): #check the bending moment of all members
    
    FirstOrderResponse, SecondOrderResponseT = setup_model

    mf = FirstOrderResponse.MemberForceLocal(1, All = True)
    BendingMomentallT =[]
    for i in range(0,9):
        MemNo = i+1
        BendingMomentallT = BendingMomentallT + FirstOrderResponse.MemberBMD(MemNo, mf[i])[:-1]
    BendingMomentallT = BendingMomentallT + FirstOrderResponse.MemberBMD(MemNo, mf[9])

    BendingMomentallR = [
     0.00, 466.20, 932.40, 1398.60, 1864.79, 2330.99, 2797.19, 3263.39, 3729.59, 4195.79,
     4661.98, 5128.18, 5594.38, 6060.58, 6526.78, 6992.98, 7459.17, 7925.37, 8391.57, 8857.77,
     9323.97, 9790.17, 10256.37, 10722.56, 11188.76, 11654.96, 12121.16, 12587.36, 13053.56, 13519.75,
     13985.95, 14452.15, 14918.35, 15384.55, 15850.75, 16316.95, 16783.14, 17249.34, 17715.54, 18181.74,
     18647.94, 19114.13, 19580.33, 20046.53, 20512.73, 20978.93, 21445.13, 21911.33, 22377.53, 22843.72,
     23309.92, 23776.12, 24242.32, 24708.52, 25174.71, 25640.91, 26107.11, 26573.31, 27039.51, 27505.71,
     27971.91, 28438.10, 28904.30, 29370.50, 29836.70, 30302.90, 30769.10, 31235.29, 31701.49, 32167.69,
     32633.89, 33100.09, 33566.29, 34032.49, 34498.68, 34964.88, 35431.08, 35897.28, 36363.48, 36829.68,
     37295.88, 37762.07, 38228.27, 38694.47, 39160.67, 39626.87, 40093.07, 40559.26, 41025.46, 41491.66,
     41957.86, 42424.06, 42890.25, 43356.46, 43822.65, 44288.85, 44755.05, 45221.25, 45687.45, 46153.64,
     46619.84
]
    # Reference from sofistik

    for i in range(100):
        assert np.allclose(BendingMomentallR[i], BendingMomentallT[i], atol=100), "BM is wrong."
    # SInce its a indeterminant structure, Moment at corner is property of stiffness, hence a tolereance of 100 for 46000 is acceptable
    #assert np.allclose(BendingMomentallR, BendingMomentallT, atol=100), "First Order BM is wrong."

    #First Order displacements have been checked with app its correct, also mmeber force locals


def test_SecondOrderBendingMomentAll(setup_model): #check the bending moment of all members
    
    FirstOrderResponse, SecondOrderResponse = setup_model

    mf = SecondOrderResponse.MemberForceLocal(1, All = True)
    BendingMomentallT =[]
    for i in range(0,9):
        MemNo = i+1
        BendingMomentallT = BendingMomentallT + SecondOrderResponse.MemberBMD(MemNo, mf[i])[:-1]
    BendingMomentallT = BendingMomentallT + SecondOrderResponse.MemberBMD(MemNo, mf[9])

    BendingMomentallR = [
    0.00, 2816.88, 5631.08, 8439.92, 11240.72, 14030.82, 16807.55, 19568.29, 22310.39, 25031.25,
    27728.28, 30398.91, 33040.59, 35650.81, 38227.09, 40766.97, 43268.03, 45727.90, 48144.23, 50514.72,
    52837.11, 55109.19, 57328.81, 59493.84, 61602.22, 63651.95, 65641.08, 67567.71, 69430.01, 71226.20,
    72954.58, 74613.49, 76201.38, 77716.70, 79158.03, 80524.00, 81813.30, 83024.71, 84157.07, 85209.31,
    86180.42, 87069.49, 87875.66, 88598.16, 89236.32, 89789.52, 90257.23, 90639.02, 90934.52, 91143.45,
    91265.60, 91300.88, 91249.22, 91110.70, 90885.45, 90573.67, 90175.67, 89691.83, 89122.60, 88468.53,
    87730.25, 86908.45, 86003.92, 85017.52, 83950.20, 82802.97, 81576.91, 80273.21, 78893.09, 77437.89,
    75908.98, 74307.82, 72635.93, 70894.90, 69086.41, 67212.16, 65273.95, 63273.61, 61213.05, 59094.24,
    56919.20, 54689.98, 52408.75, 50077.62, 47698.88, 45274.72, 42807.49, 40299.55, 37753.26, 35171.08,
    32555.42, 29908.80, 27233.76, 24532.80, 21808.55, 19063.54, 16300.41, 13521.82, 10730.37, 7928.77,
    5119.63
    ]
    # Reference from sofistik

    #for i in range(100):
    #    assert np.allclose(BendingMomentallR[i], BendingMomentallT[i], atol=2500), "BM is wrong."
    # Since its a indeterminant structure, Moment at corner is property of stiffness, 
    # Minimum amount of Finite element used,
    # hence a tolereance of 2500 for 91265 is acceptable. Pattern following checked
    assert np.allclose(BendingMomentallR, BendingMomentallT, atol=2500), "Second Order BM is wrong."







