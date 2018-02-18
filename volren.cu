#include "volren.cuh"

bool QuadNodeD_insideQuad(const QuadNodeD &q, float x, float y)
{
  return x >= q.Ax && x < q.Bx && y >= q.Ay && y < q.By;
}

bool QuadNodeD_insideTriangle(const QuadNodeD &q, float x, float y, float &alpha, float &beta, float &gamma) 
{
  alpha = ((q.y1 - q.y2)*(x - q.x2) + (q.x2 - q.x1)*(y - q.y2)) /
          ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  beta = ((q.y2 - q.y0)*(x - q.x2) + (q.x0 - q.x2)*(y - q.y2)) /
         ((q.y1 - q.y2)*(q.x0 - q.x2) + (q.x2 - q.x1)*(q.y0 - q.y2));
  gamma = 1.0 - alpha - beta;
  // fprintf(stderr, "barycentric=%f, %f, %f\n", alpha, beta, gamma);
  return alpha >= 0 && beta >= 0 && gamma >= 0;
}

int QuadNodeD_locatePoint_recursive(const QuadNodeD *q, const QuadNodeD *nodes, float x, float y, float &alpha, float &beta, float &gamma)
{
  if (q->triangleId >= 0) { //leaf node
    bool succ = QuadNodeD_insideTriangle(*q, x, y, alpha, beta, gamma);
    if (succ) return q->triangleId;
  } else if (QuadNodeD_insideQuad(*q, x, y)) {
    for (int j=0; j<4; j++) {
      if (q->childrenIds[j] > 0) {
        int result = QuadNodeD_locatePoint_recursive(&nodes[q->childrenIds[j]], nodes, x, y, alpha, beta, gamma);
        if (result >= 0) return result;
      }
    }
  }
  return -1;
}

int QuadNodeD_locatePoint(QuadNodeD *nodes, float x, float y, float &alpha, float &beta, float &gamma)
{
  // float alpha, beta, gamma;
  static const int maxStackSize = 64;
  int stack[maxStackSize];
  int stackPos = 0;
  stack[stackPos++] = 0; // push root

  while (stackPos > 0) {
    const int i = stack[--stackPos]; // pop
    const QuadNodeD &q = nodes[i];

    // fprintf(stderr, "D_checking node %d, %f, %f, %f, %f\n", i, q.Ax, q.Ay, q.Bx, q.By);
    // fprintf(stderr, "D_checking node %d\n", i);

    if (q.triangleId >= 0) { // leaf node
      bool succ = QuadNodeD_insideTriangle(q, x, y, alpha, beta, gamma);
      if (succ) return q.triangleId;
    } else if (QuadNodeD_insideQuad(q, x, y)) { // non-leaf node
      for (int j=0; j<4; j++) {
        if (q.childrenIds[j] > 0)
          stack[stackPos++] = q.childrenIds[j];
      }
    }
  }
  return -1;
}

float QuadNodeD_sample(QuadNodeD *nodes, float x, float y, float *scalar) {
  float alpha, beta, gamma;
  int i = QuadNodeD_locatePoint(nodes, x, y, alpha, beta, gamma);
  const QuadNodeD &q = nodes[i];

  return alpha * scalar[q.i0] 
    + beta * scalar[q.i1]
    + gamma * scalar[q.i2];
}
