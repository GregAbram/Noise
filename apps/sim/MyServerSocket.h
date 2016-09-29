#include <sys/select.h>
#include <vtkObjectFactory.h>
#include <vtkServerSocket.h>

class MyServerSocket : public vtkServerSocket
{
public:
	static MyServerSocket *New();
	bool ConnectionWaiting();
};

vtkStandardNewMacro(MyServerSocket);

bool 
MyServerSocket::ConnectionWaiting()
{
	int s = GetSocketDescriptor();

	fd_set rset;
	FD_ZERO(&rset);
	FD_SET(s, &rset);

	timeval tval;
	tval.tv_sec = 0;
	tval.tv_usec = 0;

	return 1 == select(s+1, &rset, NULL, NULL, &tval);
}

