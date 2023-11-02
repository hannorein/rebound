/* 
 * tiny.c - a minimal HTTP server
 *          Dave O'Hallaron, Carnegie Mellon
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <netdb.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include "rebound.h"

#define BUFSIZE 1024

extern char **environ; /* the environment */

/*
 * cerror - returns an error message to the client
 */
void cerror(FILE *stream, char *cause, char *errno, 
        char *shortmsg, char *longmsg) {
    fprintf(stream, "HTTP/1.1 %s %s\n", errno, shortmsg);
    fprintf(stream, "Content-type: text/html\n");
    fprintf(stream, "\n");
    fprintf(stream, "<html><title>Tiny Error</title>");
    fprintf(stream, "<body bgcolor=""ffffff"">\n");
    fprintf(stream, "%s: %s\n", errno, shortmsg);
    fprintf(stream, "<p>%s: %s\n", longmsg, cause);
    fprintf(stream, "<hr><em>The Tiny Web server</em>\n");
}

void* start_server(void* args){
    pthread_setcancelstate(PTHREAD_CANCEL_ENABLE, NULL);
    pthread_setcanceltype(PTHREAD_CANCEL_DEFERRED, NULL);
    struct reb_simulation* r = (struct reb_simulation*)args;

    /* variables for connection management */
    int parentfd;          /* parent socket */
    int childfd;           /* child socket */
    int portno;            /* port to listen on */
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

    /* check command line args */
    portno = 1234;

    /* open socket descriptor */
    parentfd = socket(AF_INET, SOCK_STREAM, 0);
    if (parentfd < 0) 
        reb_exit("ERROR opening socket");

    /* allows us to restart server immediately */
    optval = 1;
    setsockopt(parentfd, SOL_SOCKET, SO_REUSEADDR, 
            (const void *)&optval , sizeof(int));

    /* bind port to socket */
    bzero((char *) &serveraddr, sizeof(serveraddr));
    serveraddr.sin_family = AF_INET;
    serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
    serveraddr.sin_port = htons((unsigned short)portno);
    if (bind(parentfd, (struct sockaddr *) &serveraddr, 
                sizeof(serveraddr)) < 0) 
        reb_exit("ERROR on binding");

    /* get us ready to accept connection requests */
    if (listen(parentfd, 5) < 0) /* allow 5 requests to queue up */ 
        reb_exit("ERROR on listen");

    printf("Listening on port %d...\n",portno);
    /* 
     * main loop: wait for a connection request, parse HTTP,
     * serve requested content, close connection.
     */
    clientlen = sizeof(clientaddr);
    while (1) {
        /* wait for a connection request */
        printf("Pre accept\n");
        childfd = accept(parentfd, (struct sockaddr *) &clientaddr, &clientlen);
        printf("Post accept\n");
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
            cerror(stream, method, "501", "Not Implemented", "Tiny did not get request");
            fclose(stream);
            close(childfd);
            continue;
        }

        sscanf(buf, "%s %s %s\n", method, uri, version);

        /* tiny only supports the GET method */
        if (strcasecmp(method, "GET")) {
            cerror(stream, method, "501", "Not Implemented", "Tiny does not implement this method");
            fclose(stream);
            close(childfd);
            continue;
        }
            
        printf("URI: %s\n",uri);

        /* read (and ignore) the HTTP headers */
        fgets(buf, BUFSIZE, stream);
        //printf("%s", buf);
        while(strcmp(buf, "\r\n")) {
            fgets(buf, BUFSIZE, stream);
            //printf("%s", buf);
        }

        fprintf(stream, "HTTP/1.1 200 OK\n");
        fprintf(stream, "Server: Tiny Web Server\n");
        if (!strcasecmp(uri, "/simulation")) {
            char* bufp = NULL;
            size_t sizep;
            reb_simulation_save_to_stream(r, &bufp,&sizep);
            //fprintf(stream, "Content-length: %d\n", (int)sizep);
            //fprintf(stream, "Content-type: application/octet-stream\n");
            fprintf(stream, "Content-type: text/html\n");
            fprintf(stream, "Access-Control-Allow-Origin: *\n");
            fprintf(stream, "Cross-Origin-Opener-Policy: cross-origin\n");
            fprintf(stream, "\r\n"); 
            fflush(stream);
            fwrite(bufp, 1, sizep, stream);
            free(bufp);
        }else{
            fprintf(stream, "Content-type: text/html\n");
            fprintf(stream, "\r\n"); 
            fflush(stream);
            fprintf(stream, "<h1>Not sure what to do!</h1>\n");
        }

        /* clean up */
        fclose(stream);
        close(childfd);

    }
    printf("Server shutting down...\n");
}
