#ifndef enum_h
#define enum_h

namespace OF {

namespace Geometry {

enum  Sign
{
    NEGATIVE = -1, ZERO = 0, POSITIVE = 1,

    // Orientation constants:
    RIGHT_TURN = -1, LEFT_TURN = 1,

    CLOCKWISE = -1, COUNTERCLOCKWISE = 1,

    COLLINEAR = 0, COPLANAR = 0, DEGENERATE = 0,

    // Oriented_side constants:
    ON_NEGATIVE_SIDE = -1, ON_ORIENTED_BOUNDARY = 0, ON_POSITIVE_SIDE = 1,

    // Comparison_result constants:
    SMALLER = -1, EQUAL = 0, LARGER = 1
};

enum  Bounded_side
      {
        ON_UNBOUNDED_SIDE = -1,
        ON_BOUNDARY,
        ON_BOUNDED_SIDE
      };

enum  Angle
      {
          OBTUSE = -1,
          RIGHT,
          ACUTE
      };

} // end of namespace Geometry
} // end of namespace OF
#endif // end of enum_h
