#include "concurrentqueue.h"

int main(int argc, char **argv) 
{
  moodycamel::ConcurrentQueue<int> q;
  q.enqueue(25);

  int item;
  bool found = q.try_dequeue(item);
  assert(found && item == 25);

  exit(0);
}
