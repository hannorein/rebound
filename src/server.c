/**
 * @file    server.c
 * @brief   Opens a webserver to allow for platform independent visualization.
 * @author  Hanno Rein <hanno@hanno-rein.de>
 * @details These functions provide real time visualizations
 * using OpenGL. Part of the code is by Dave O'Hallaron, Carnegie Mellon (tiny.c).
 * 
 * @section LICENSE
 * Copyright (c) 2023 Hanno Rein, Dave O'Hallaron, Carnegie Mellon
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#ifdef _MSC_VER 
//not #if defined(_WIN32) || defined(_WIN64) because we have strncasecmp in mingw
#define strncasecmp _strnicmp
#define strcasecmp _stricmp
#endif
#ifdef _WIN32
#include <WinSock2.h>
#include <WS2tcpip.h>
#include <tchar.h>

// Buffer 
#define Buffer 1500
// A Pragma comment. Loading Library 
#pragma comment(lib, "ws2_32.lib")
#else // _WIN32
#include <unistd.h>
#include <netdb.h>
#include <sys/socket.h>
#include <sys/mman.h>
#include <netinet/in.h>
#endif // _WIN32
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include "rebound.h"

#define BUFSIZE 1024

static void reb_server_cerror(FILE *stream, char *cause, char *errno, char *shortmsg, char *longmsg) {
    fprintf(stream, "HTTP/1.1 %s %s\n", errno, shortmsg);
    fprintf(stream, "Content-type: text/html\n");
    fprintf(stream, "\n");
    fprintf(stream, "<html><title>REBOUND Webserver Error</title>");
    fprintf(stream, "<body bgcolor=""ffffff"">\n");
    fprintf(stream, "%s: %s\n", errno, shortmsg);
    fprintf(stream, "<p>%s: %s\n", longmsg, cause);
    fprintf(stream, "<hr><em>REBOUND Webserver</em>\n");
}

void* reb_server_start(void* args){
#ifndef _WIN32
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);

    struct reb_server_data* data = (struct reb_server_data*)args;
    struct reb_simulation* r = data->r;

    /* variables for connection management */
    int parentfd;          /* parent socket */
    int childfd;           /* child socket */
    unsigned int clientlen;         /* byte size of client's address */
    struct hostent *hostp; /* client host info */
    int optval;            /* flag value for setsockopt */
    struct sockaddr_in serveraddr; /* server's addr */
    struct sockaddr_in clientaddr; /* client addr */

    /* variables for connection I/O */
    FILE *stream;          /* stream version of childfd */
    char buf[BUFSIZE];     /* message buffer */
    char method[BUFSIZE];  /* request method */
    char uri[BUFSIZE];     /* request uri */
    char version[BUFSIZE]; /* request method */

    /* open socket descriptor */
    parentfd = socket(AF_INET, SOCK_STREAM, 0);
    if (parentfd < 0)
        reb_exit("ERROR opening socket");

    /* allows us to restart server immediately */
    optval = 1;
    setsockopt(parentfd, SOL_SOCKET, SO_REUSEADDR,
            (const void *)&optval , sizeof(int));

    /* bind port to socket */
    memset((char *) &serveraddr, 0, sizeof(serveraddr));
    serveraddr.sin_family = AF_INET;
    serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
    serveraddr.sin_port = htons((unsigned short)data->port);
    if (bind(parentfd, (struct sockaddr *) &serveraddr, sizeof(serveraddr)) < 0){
        printf("Error opening binding port %d. Port might be in use.\n", data->port);
        return PTHREAD_CANCELED;
    }

    /* get us ready to accept connection requests */
    if (listen(parentfd, 5) < 0) /* allow 5 requests to queue up */
        reb_exit("ERROR on listen");

    printf("REBOUND Webserver listening on http://localhost:%d ...\n",data->port);
    /*
     * main loop: wait for a connection request, parse HTTP,
     * serve requested content, close connection.
     */
    clientlen = sizeof(clientaddr);
    while (1) {
        /* wait for a connection request */
        childfd = accept(parentfd, (struct sockaddr *) &clientaddr, &clientlen);
        if (childfd < 0)
            reb_exit("ERROR on accept");

        /* determine who sent the message */
        hostp = gethostbyaddr((const char *)&clientaddr.sin_addr.s_addr,
                sizeof(clientaddr.sin_addr.s_addr), AF_INET);
        if (hostp == NULL)
            reb_exit("ERROR on gethostbyaddr");

        /* open the child socket descriptor as a stream */
        if ((stream = fdopen(childfd, "r+")) == NULL)
            reb_exit("ERROR on fdopen");

        /* get the HTTP request line */
        char* request = fgets(buf, BUFSIZE, stream);
        if (!request){
            reb_server_cerror(stream, method, "501", "Not Implemented", "REBOUND Webserver did not get request");
            fclose(stream);
            close(childfd);
            continue;
        }

        sscanf(buf, "%s %s %s\n", method, uri, version);

        /* only support the GET method */
        if (strcasecmp(method, "GET")) {
            reb_server_cerror(stream, method, "501", "Not Implemented", "REBOUND Webserver does not implement this method");
            fclose(stream);
            close(childfd);
            continue;
        }
           

        /* read (and ignore) the HTTP headers */
        fgets(buf, BUFSIZE, stream);
        while(strcmp(buf, "\r\n")) {
            fgets(buf, BUFSIZE, stream);
        }

        fprintf(stream, "HTTP/1.1 200 OK\n");
        fprintf(stream, "Server: REBOUND Webserver\n");
        //fprintf(stream, "Access-Control-Allow-Origin: *\n");
        //fprintf(stream, "Cross-Origin-Opener-Policy: cross-origin\n");
        fprintf(stream, "Content-type: text/html\n"); // Always using the same content type, even for binary data.
        fprintf(stream, "\r\n");
        if (!strcasecmp(uri, "/simulation")) {
            char* bufp = NULL;
            size_t sizep;
            data->need_copy = 1;
            pthread_mutex_lock(&(data->mutex));
            reb_simulation_save_to_stream(r, &bufp,&sizep);
            data->need_copy = 0;
            pthread_mutex_unlock(&(data->mutex));
            fflush(stream);
            fwrite(bufp, 1, sizep, stream);
            free(bufp);
        }else if (!strncasecmp(uri, "/keyboard/",10)) {
            int key = 0;
            sscanf(uri, "/keyboard/%d", &key);
            switch (key){
                case 'Q':
                    data->r->status = REB_STATUS_USER;
                    fprintf(stream, "ok.\n");
                    break;
                case ' ':
                    if (data->r->status == REB_STATUS_PAUSED){
                        printf("Resume.\n");
                        data->r->status = REB_STATUS_RUNNING;
                    }else{
                        printf("Pause.\n");
                        data->r->status = REB_STATUS_PAUSED;
                    }
                    fprintf(stream, "ok.\n");
                    break;
                default:
                    fprintf(stream, "Unknown key received: %d\n",key);
                    break;
            }
            fflush(stream);
        }else if (!strcasecmp(uri, "/") || !strcasecmp(uri, "/index.html") || !strcasecmp(uri, "/rebound.html")) {
            struct stat sbuf;
            if (stat("rebound.html", &sbuf) < 0) {
                printf("Unable to find rebound.html\n");
                fflush(stream);
                fprintf(stream, "<h1>Unable to find rebound.html.</h1>\n");
                fprintf(stream, "The server was unable to find rebound.html.\n");
                fprintf(stream, "Download rebound.html from github and place it in the same directory as your rebound executable.\n");
            }else{
                fflush(stream);
#ifndef _WIN32
                int fd = open("rebound.html", O_RDONLY);
                void* p = mmap(0, sbuf.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
                fwrite(p, 1, sbuf.st_size, stream);
                munmap(p, sbuf.st_size);
#else //_WIN32
            fprintf(stream, "Not yet implemented.\n");
#endif //_WIN32
            }
        }else if (!strcasecmp(uri, "/favicon.ico")) {
            fflush(stream);
            fprintf(stream, "favicon.ico not found.\n");
        }else{
            printf("Not sure what to do with URI: %s\n",uri);
            fflush(stream);
            fprintf(stream, "<h1>Not sure what to do!</h1>\n");
        }

        /* clean up */
        fclose(stream);
        close(childfd);

    }
    printf("Server shutting down...\n");
    return PTHREAD_CANCELED;

#else // _WIN32


    WSADATA wsa;
    struct sockaddr_in server;
    SOCKET s;
    SOCKET clientS;
    int iResult;
    const char message[Buffer] =
        "HTTP/1.1 200 OK\n"
        "Server: REBOUND Webserver\n"
        "Content-type: text/html\n"
        "\r\n"
        "Hello\n"
        ;
    char reply[Buffer] = { 0 };

    // Simple. Start WSA(Windows Sockets API). If the return answer is not 0. It means error so therefore,
    if (WSAStartup(MAKEWORD(2, 2), &wsa) != 0) {
        //std::cerr << "WinsocK Failed" << std::endl; // Print out this on the Screen.
        getchar(); // Pause Console
        exit(1); // Exit Program, Note that in Computing, 1 = True, 0 = False. So 0 or 1 can also be entered as, exit(1);
    }

    printf( "Windows Socket API Started\n");//

    // Create a Network Socket
    s = socket(AF_INET, SOCK_STREAM, NULL);
    // If the Socket is Invalid or a Socket Error Occurs
    if (s == SOCKET_ERROR)
    {
        //std::cerr << "Socket Error : " << WSAGetLastError() << std::endl;
        getchar(); // Pause Console
        exit(1); // Exit Program, Note that in Computing, 1 = True, 0 = False. So 0 or 1 can also be entered as, exit(1);
    }
    else if (s == INVALID_SOCKET) {
        //std::cerr << "Socket Error : " << WSAGetLastError() << std::endl;
        getchar(); // Pause Console
        exit(1); // Exit Program, Note that in Computing, 1 = True, 0 = False. So 0 or 1 can also be entered as, exit(1);
    }

    printf("Socket Created\n");

    server.sin_family = AF_INET; // Using AF_INET Address Family.
    server.sin_port = htons(1234); // Defining PORT
    InetPton(AF_INET, _T("0.0.0.0"), &server.sin_addr); // Defining The Network Address to Run the Server on

    iResult = bind(s, (struct sockaddr*)&server, sizeof(server)); // binding the Host Address and Port Number
    if (iResult == SOCKET_ERROR) // If Bind gives Error
    {
        //std::cerr << "Bind Error " << WSAGetLastError() << std::endl;
        getchar(); // Pause Console
        exit(1); // Exit Program, Note that in Computing, 1 = True, 0 = False. So 0 or 1 can also be entered as, exit(1);
    }

    printf("listening\n");
    //std::cout << "Listening on : 0.0.0.0:80" << std::endl; // Tell the User we Started Listening.
    iResult = listen(s, AF_INET); // Then Actually Start Listening for incoming Connections.

    /*
       The Program will start to listen for incoming connections and will do so until
       Someone Connects.
       */

    while(1){
        clientS = accept(s, NULL, NULL); // Accept a Connection on a new Socket made for the Client.
        if (clientS == SOCKET_ERROR) { // if Accepting Connection is a Error
                                       //std::cerr << "Accept FAiled!" << WSAGetLastError() << std::endl;
            getchar(); // Pause Console
            exit(1); // Exit Program, Note that in Computing, 1 = True, 0 = False. So 0 or 1 can also be entered as, exit(1);
        }
        else {
            printf( "A Client Connected. Sending a Message and closing Connection\n");
            send(clientS, message, strlen(message), NULL); // Send Client a Message
        }

        recv(clientS, reply, sizeof(reply), NULL); // Just in case if the Client sends something, We Receive it.
        printf("closing1\n");

        closesocket(clientS); // close the Client Socket now that our Work is Complete.
        printf("closing2\n");
    }
    WSACleanup(); // Clean Windows Socket API.
    printf("closing3\n");

    return NULL;
#endif // _WIN32
}
