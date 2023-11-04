/*
 * tiny.c - a minimal HTTP server
 *          Dave O'Hallaron, Carnegie Mellon
 */

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <netdb.h>
#include <fcntl.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/mman.h>
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
    bzero((char *) &serveraddr, sizeof(serveraddr));
    serveraddr.sin_family = AF_INET;
    serveraddr.sin_addr.s_addr = htonl(INADDR_ANY);
    serveraddr.sin_port = htons((unsigned short)data->port);
    if (bind(parentfd, (struct sockaddr *) &serveraddr,
                sizeof(serveraddr)) < 0)
        reb_exit("ERROR on binding");

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
           

        /* read (and ignore) the HTTP headers */
        fgets(buf, BUFSIZE, stream);
        while(strcmp(buf, "\r\n")) {
            fgets(buf, BUFSIZE, stream);
        }

        fprintf(stream, "HTTP/1.1 200 OK\n");
        fprintf(stream, "Server: Tiny Web Server\n");
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
                int fd = open("rebound.html", O_RDONLY);
                void* p = mmap(0, sbuf.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
                fwrite(p, 1, sbuf.st_size, stream);
                munmap(p, sbuf.st_size);
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
}
