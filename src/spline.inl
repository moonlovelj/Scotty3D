// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
      const T& position0,
      const T& position1,
      const T& tangent0,
      const T& tangent1,
      double normalizedTime,
      int derivative )
{
   // TODO IMPLEMENT ME (TASK 1A)
  if (derivative == 2) {
    return (12 * normalizedTime - 6) * position0 +
           (6 * normalizedTime - 4) * tangent0 +
           (-12 * normalizedTime + 6) * position1 +
           (6 * normalizedTime - 2) * tangent1;
  }
  else if (derivative == 1) {
    double t1 = normalizedTime;
    double t2 = normalizedTime * normalizedTime;
    return (6 * t2 - 6 * t1) * position0 + (3 * t2 - 4 * t1 + 1) * tangent0 +
           (-6 * t2 + 6 * t1) * position1 + (3 * t2 - 2 * t1) * tangent1;
  }
  else {
    double t1 = normalizedTime;
    double t2 = normalizedTime * normalizedTime;
    double t3 = normalizedTime * normalizedTime * normalizedTime;
    return (2 * t3 - 3 * t2 + 1) * position0 + (t3 - 2 * t2 + t1) * tangent0 +
           (-2 * t3 + 3 * t2) * position1 + (t3 - t2) * tangent1;
  }
}
            
// Returns a state interpolated between the values directly before and after the given time.
template <class T>
inline T Spline<T>::evaluate( double time, int derivative )
{
   // TODO IMPLEMENT ME (TASK 1B)
  if (knots.size() < 1) {
    return T();
  } else if (knots.size() == 1) {
    return (derivative == 0 ? knots.begin()->second : 0);
  } else {
    if (time <= knots.begin()->first) {
      return (derivative == 0 ? knots.begin()->second : T());
    } else if (time >= knots.rbegin()->first) {
      return (derivative == 0 ? knots.rbegin()->second : T());
    } else {
      std::map<double, T>::iterator t2_iter = knots.upper_bound(time);
      std::map<double, T>::iterator t1_iter = t2_iter;
      t1_iter--;
      double t1 = t1_iter->first;
      double t2 = t2_iter->first;
      T position1 = t1_iter->second;
      T position2 = t2_iter->second;
      double t0, t3;
      T position0, position3;
      if (t1_iter == knots.begin()) {
        t0 = t1 - (t2 - t1);
        position0 = evaluate(t0, derivative);
      } else {
        t0 = (--t1_iter)->first;
        position0 = t1_iter->second;
      }
      std::map<double, T>::reverse_iterator t2_reverse_iter(t2_iter);
      if (t2_reverse_iter == knots.rbegin()) {
        t3 = t2 + (t2 - t1);
        position3 = evaluate(t3, derivative);
      } else {
        t3 = (++t2_iter)->first;
        position3 = t2_iter->second;
      }
      T m1 = (position2 - position0) / (t2 - t0);
      T m2 = (position3 - position1) / (t3 - t1);
      double normalized_time = (time - t1) / (t2 - t1);
      return cubicSplineUnitInterval(position1, position2, m1, m2, normalized_time, derivative);
    }
  }
}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance )
{
   // Empty maps have no knots.
   if( knots.size() < 1 )
   {
      return false;
   }

   // Look up the first element > or = to time.
   typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
   typename std::map<double, T>::iterator t1_iter;
   t1_iter = t2_iter;
   t1_iter--;

   if( t2_iter == knots.end() )
   {
      t2_iter = t1_iter;
   }

   // Handle tolerance bounds,
   // because we are working with floating point numbers.
   double t1 = (*t1_iter).first;
   double t2 = (*t2_iter).first;

   double d1 = fabs(t1 - time);
   double d2 = fabs(t2 - time);


   if(d1 < tolerance && d1 < d2)
   {
      knots.erase(t1_iter);
      return true;
   }

   if(d2 < tolerance && d2 < d1)
   {
      knots.erase(t2_iter);
      return t2;
   }

   return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue( double time, T value )
{
   knots[ time ] = value;
}

template <class T>
inline T Spline<T>::operator()( double time )
{
   return evaluate( time );
}
