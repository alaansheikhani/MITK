/*=========================================================================

Program:   openCherry Platform
Language:  C++
Date:      $Date$
Version:   $Revision$

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/

#ifndef CHERRYMESSAGE_H_
#define CHERRYMESSAGE_H_

#include <vector>
#include <Poco/Mutex.h>

#define NewMessageMacro(msgName)                                        \
  public:                                                               \
    void AddmsgNameListener(const MessageAbstractDelegate& delegate)    \
    { m_msgName += delegate; }                                          \
    void RemovemsgNameListener(const MessageAbstractDelegate& delegate) \
    { m_msgName -= delegate; }                                          \
  private: Message m_msgName;                                           \

#define NewMessage1Macro(msgName, type1)                                          \
  public:                                                                         \
    void AddmsgNameListener(const MessageAbstractDelegate1< type1 >& delegate)    \
    { m_msgName += delegate; }                                                    \
    void RemovemsgNameListener(const MessageAbstractDelegate1< type1 >& delegate) \
    { m_msgName -= delegate; }                                                    \
  private: Message1< type1 > m_msgName;                                           \


namespace cherry {

class MessageAbstractDelegate
{
  public:

    virtual ~MessageAbstractDelegate()
    {
    }

    virtual void Execute() = 0;
    virtual bool operator==(const MessageAbstractDelegate* cmd) = 0;
    virtual MessageAbstractDelegate* Clone() const = 0;
};


template <typename T, typename A = void>
class MessageAbstractDelegate1
{
  public:

    virtual ~MessageAbstractDelegate1()
    {
    }

    virtual A Execute(T t) = 0;
    virtual bool operator==(const MessageAbstractDelegate1* cmd) = 0;
    virtual MessageAbstractDelegate1* Clone() const = 0;
};

template <typename T, typename U, typename A = void>
class MessageAbstractDelegate2
{
  public:

    virtual ~MessageAbstractDelegate2()
    {
    }

    virtual A Execute(T t, U u) const = 0;
    virtual bool operator==(const MessageAbstractDelegate2* cmd) = 0;
    virtual MessageAbstractDelegate2* Clone() const = 0;
};

template <typename T, typename U, typename V, typename A = void>
class MessageAbstractDelegate3
{
  public:

    virtual ~MessageAbstractDelegate3()
    {
    }

    virtual A Execute(T t, U u, V v) const = 0;
    virtual bool operator==(const MessageAbstractDelegate3* cmd) = 0;
    virtual MessageAbstractDelegate3* Clone() const = 0;
};

template <typename T, typename U, typename V, typename W, typename A = void>
class MessageAbstractDelegate4
{
  public:

    virtual ~MessageAbstractDelegate4()
    {
    }

    virtual A Execute(T t, U u, V v, W w) const = 0;
    virtual bool operator==(const MessageAbstractDelegate4* cmd) = 0;
    virtual MessageAbstractDelegate4* Clone() const = 0;
};

template <class R>
class MessageDelegate : public MessageAbstractDelegate
{
  public:

    // constructor - takes pointer to an object and pointer to a member and stores
    // them in two private variables
    MessageDelegate(R* object, void(R::*memberFunctionPointer)())
    :m_Object(object),
    m_MemberFunctionPointer(memberFunctionPointer)
    {
    }

    virtual ~MessageDelegate()
    {
    }

    // override function "Call"
    virtual void Execute()
    {
      return (m_Object->*m_MemberFunctionPointer)();    // execute member function
    }

    bool operator==(const MessageAbstractDelegate* c)
    {
      const MessageDelegate<R>* cmd = dynamic_cast<const MessageDelegate<R>* >(c);
      if (!cmd) return false;

      if ((void*)this->m_Object != (void*)cmd->m_Object) return false;
      if (this->m_MemberFunctionPointer != cmd->m_MemberFunctionPointer) return false;
      return true;
    }

    MessageAbstractDelegate* Clone() const
    {
      return new MessageDelegate(m_Object, m_MemberFunctionPointer);
    }

  private:
    R* m_Object;                            // pointer to object
    void (R::*m_MemberFunctionPointer)();   // pointer to member function
};


template <class R, typename T, typename A = void>
class MessageDelegate1 : public MessageAbstractDelegate1<T,A>
{
  public:

    // constructor - takes pointer to an object and pointer to a member and stores
    // them in two private variables
    MessageDelegate1(R* object, A(R::*memberFunctionPointer)(T))
    :m_Object(object),
    m_MemberFunctionPointer(memberFunctionPointer)
    {
    }

    virtual ~MessageDelegate1()
    {
    }

    // override function "Call"
    virtual A Execute(T t)
    {
      return (m_Object->*m_MemberFunctionPointer)(t);    // execute member function
    }

    bool operator==(const MessageAbstractDelegate1<T,A>* c)
    {
      const MessageDelegate1<R,T,A>* cmd = dynamic_cast<const MessageDelegate1<R,T,A>* >(c);
      if (!cmd) return false;

      if ((void*)this->m_Object != (void*)cmd->m_Object) return false;
      if (this->m_MemberFunctionPointer != cmd->m_MemberFunctionPointer) return false;
      return true;
    }

    MessageAbstractDelegate1<T,A>* Clone() const
    {
      return new MessageDelegate1(m_Object, m_MemberFunctionPointer);
    }

  private:
    R* m_Object;                            // pointer to object
    A (R::*m_MemberFunctionPointer)(T);   // pointer to member function
};


template <class R, typename T, typename U, typename A = void>
class MessageDelegate2 : public MessageAbstractDelegate2<T,U,A>
{
  public:

    // constructor - takes pointer to an object and pointer to a member and stores
    // them in two private variables
    MessageDelegate2(R* object, A(R::*memberFunctionPointer)(T, U))
    :m_Object(object),
    m_MemberFunctionPointer(memberFunctionPointer)
    {
    }

    virtual ~MessageDelegate2()
    {
    }

    // override function "Call"
    virtual A Execute(T t, U u) const
    {
      return (m_Object->*m_MemberFunctionPointer)(t,u);             // execute member function
    }

    bool operator==(const MessageAbstractDelegate2<T,U,A>* c)
    {
      const MessageDelegate2<R,T,U,A>* cmd = dynamic_cast<const MessageDelegate2<R,T,U,A>* >(c);
      if (!cmd) return false;

      if ((void*)this->m_Object != (void*)cmd->m_Object) return false;
      if (this->m_MemberFunctionPointer != cmd->m_MemberFunctionPointer) return false;
      return true;
    }

    MessageAbstractDelegate2<T,U,A>* Clone() const
    {
      return new MessageDelegate2(m_Object, m_MemberFunctionPointer);
    }

  private:
    R* m_Object;                            // pointer to object
    A (R::*m_MemberFunctionPointer)(T, U);   // pointer to member function
};

template <class R, typename T, typename U, typename V, typename A = void>
class MessageDelegate3 : public MessageAbstractDelegate3<T,U,V,A>
{
  public:

    // constructor - takes pointer to an object and pointer to a member and stores
    // them in two private variables
    MessageDelegate3(R* object, A(R::*memberFunctionPointer)(T, U, V))
    :m_Object(object),
    m_MemberFunctionPointer(memberFunctionPointer)
    {
    }

    virtual ~MessageDelegate3()
    {
    }

    // override function "Call"
    virtual A Execute(T t, U u, V v) const
    {
      return (m_Object->*m_MemberFunctionPointer)(t,u,v);             // execute member function
    }

    bool operator==(const MessageAbstractDelegate3<T,U,V,A>* c)
    {
      const MessageDelegate3<R,T,U,V,A>* cmd = dynamic_cast<const MessageDelegate3<R,T,U,V,A>* >(c);
      if (!cmd) return false;

      if ((void*)this->m_Object != (void*)cmd->m_Object) return false;
      if (this->m_MemberFunctionPointer != cmd->m_MemberFunctionPointer) return false;
      return true;
    }

    MessageAbstractDelegate3<T,U,V,A>* Clone() const
    {
      return new MessageDelegate3(m_Object, m_MemberFunctionPointer);
    }

  private:
    R* m_Object;                            // pointer to object
    A (R::*m_MemberFunctionPointer)(T, U, V);   // pointer to member function
};

template <class R, typename T, typename U, typename V, typename W, typename A = void>
class MessageDelegate4 : public MessageAbstractDelegate4<T,U,V,W,A>
{
  public:

    // constructor - takes pointer to an object and pointer to a member and stores
    // them in two private variables
    MessageDelegate4(R* object, A(R::*memberFunctionPointer)(T, U, V, W))
    :m_Object(object),
    m_MemberFunctionPointer(memberFunctionPointer)
    {
    }

    virtual ~MessageDelegate4()
    {
    }

    // override function "Call"
    virtual A Execute(T t, U u, V v, W w) const
    {
      return (m_Object->*m_MemberFunctionPointer)(t,u,v,w);             // execute member function
    }

    bool operator==(const MessageAbstractDelegate4<T,U,V,W,A>* c)
    {
      const MessageDelegate4<R,T,U,V,W,A>* cmd = dynamic_cast<const MessageDelegate4<R,T,U,V,W,A>* >(c);
      if (!cmd) return false;

      if ((void*)this->m_Object != (void*)cmd->m_Object) return false;
      if (this->m_MemberFunctionPointer != cmd->m_MemberFunctionPointer) return false;
      return true;
    }

    MessageAbstractDelegate4<T,U,V,W,A>* Clone() const
    {
      return new MessageDelegate4(m_Object, m_MemberFunctionPointer);
    }

  private:
    R* m_Object;                            // pointer to object
    A (R::*m_MemberFunctionPointer)(T, U, V, W);   // pointer to member function
};

// message without parameters (pure signals)
class Message
{
  public:

    typedef Message Self;
    typedef MessageAbstractDelegate AbstractDelegate;
    typedef std::vector<AbstractDelegate* > ListenerList;

    ~Message() {
      for (ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end(); ++iter )
      {
        delete *iter;
      }
    }

    void AddListener( const AbstractDelegate& delegate ) const
    {
      AbstractDelegate* msgCmd = delegate.Clone();

      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(msgCmd)) {
          delete msgCmd;
          return;
        }
      }
      m_Listeners.push_back(msgCmd);
    }

    void operator += ( const AbstractDelegate& delegate ) const
    {
      this->AddListener(delegate);
    }

    void RemoveListener( const AbstractDelegate& delegate ) const
    {
      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(&delegate))
        {
          delete *iter;
          m_Listeners.erase( iter );
          return;
        }
      }
    }

    void operator -= ( const AbstractDelegate& delegate) const
    {
      this->RemoveListener(delegate);
    }

    void Send()
    {
      ListenerList listeners;

      {
        Poco::FastMutex::ScopedLock lock(m_Mutex);
        listeners.assign(m_Listeners.begin(), m_Listeners.end());
      }

      for ( ListenerList::iterator iter = listeners.begin();
            iter != listeners.end();
            ++iter )
      {
        // notify each listener
        (*iter)->Execute();
      }
    }

    void operator()()
    {
      this->Send();
    }

  protected:

    mutable ListenerList m_Listeners;
    mutable Poco::FastMutex m_Mutex;

};


// message with 1 parameter and return type
template <typename T, typename A = void>
class Message1
{
  public:

    typedef Message1 Self;
    typedef MessageAbstractDelegate1<T,A> AbstractDelegate;
    typedef std::vector<AbstractDelegate*> ListenerList;

    ~Message1() {
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end(); ++iter )
      {
        delete *iter;
      }
    }

    void AddListener( const AbstractDelegate& delegate ) const
    {
      AbstractDelegate* msgCmd = delegate.Clone();

      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(msgCmd)) {
          delete msgCmd;
          return;
        }
      }
      m_Listeners.push_back(msgCmd);
    }

    void operator += ( const AbstractDelegate& delegate ) const
    {
      this->AddListener(delegate);
    }

    void RemoveListener( const AbstractDelegate& delegate ) const
    {
      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(&delegate))
        {
          delete *iter;
          m_Listeners.erase( iter );
          return;
        }
      }
    }

    void operator -= ( const AbstractDelegate& delegate) const
    {
      this->RemoveListener(delegate);
    }

    void Send(T t)
    {
      ListenerList listeners;

      {
        Poco::FastMutex::ScopedLock lock(m_Mutex);
        listeners.assign(m_Listeners.begin(), m_Listeners.end());
      }

      for ( typename ListenerList::iterator iter = listeners.begin();
            iter != listeners.end();
            ++iter )
      {
        // notify each listener
        (*iter)->Execute(t);
      }
    }

    void operator()(T t)
    {
      this->Send(t);
    }

    const ListenerList& GetListeners()
    {
      return m_Listeners;
    }

  protected:

    mutable ListenerList m_Listeners;
    mutable Poco::FastMutex m_Mutex;

};


// message with 2 parameters and return type
template <typename T, typename U, typename A = void>
class Message2
{
  public:

    typedef Message2 Self;
    typedef MessageAbstractDelegate2<T,U,A> AbstractDelegate;
    typedef std::vector<AbstractDelegate* > ListenerList;

    ~Message2() {
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end(); ++iter )
      {
        delete *iter;
      }
    }

    void AddListener( const AbstractDelegate& delegate ) const
    {
      AbstractDelegate* msgCmd = delegate.Clone();

      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(msgCmd)) {
          delete msgCmd;
          return;
        }
      }
      m_Listeners.push_back(msgCmd);
    }

    void operator += ( const AbstractDelegate& delegate ) const
    {
      this->AddListener(delegate);
    }

    void RemoveListener( const AbstractDelegate& delegate ) const
    {
      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(&delegate))
        {
          delete *iter;
          m_Listeners.erase( iter );
          return;
        }
      }
    }

    void operator -= ( const AbstractDelegate& delegate) const
    {
      this->RemoveListener(delegate);
    }

    void Send(T t, U u)
    {
      ListenerList listeners;

      {
        Poco::FastMutex::ScopedLock lock(m_Mutex);
        listeners.assign(m_Listeners.begin(), m_Listeners.end());
      }

      for ( typename ListenerList::iterator iter = listeners.begin();
            iter != listeners.end();
            ++iter )
      {
        // notify each listener
        (*iter)->Execute(t,u);
      }
    }

    void operator()(T t, U u)
    {
      this->Send(t, u);
    }

    const ListenerList& GetListeners()
    {
      return m_Listeners;
    }


  protected:

    mutable ListenerList m_Listeners;
    mutable Poco::FastMutex m_Mutex;

};

// message with 3 parameters and return type
template <typename T, typename U, typename V, typename A = void>
class Message3
{
  public:

    typedef Message3 Self;
    typedef MessageAbstractDelegate3<T,U,V,A> AbstractDelegate;
    typedef std::vector<AbstractDelegate* > ListenerList;

    ~Message3() {
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end(); ++iter )
      {
        delete *iter;
      }
    }

    void AddListener( const AbstractDelegate& delegate ) const
    {
      AbstractDelegate* msgCmd = delegate.Clone();

      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(msgCmd)) {
          delete msgCmd;
          return;
        }
      }
      m_Listeners.push_back(msgCmd);
    }

    void operator += ( const AbstractDelegate& delegate ) const
    {
      this->AddListener(delegate);
    }

    void RemoveListener(const AbstractDelegate& delegate ) const
    {
      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(&delegate))
        {
          delete *iter;
          m_Listeners.erase( iter );
          return;
        }
      }
    }

    void operator -= ( const AbstractDelegate& delegate) const
    {
      this->RemoveListener(delegate);
    }

    void Send(T t, U u, V v)
    {
      ListenerList listeners;

      {
        Poco::FastMutex::ScopedLock lock(m_Mutex);
        listeners.assign(m_Listeners.begin(), m_Listeners.end());
      }

      for ( typename ListenerList::iterator iter = listeners.begin();
            iter != listeners.end();
            ++iter )
      {
        // notify each listener
        (*iter)->Execute(t,u,v);
      }
    }

    void operator()(T t, U u, V v)
    {
      this->Send(t, u, v);
    }

    const ListenerList& GetListeners()
    {
      return m_Listeners;
    }


  protected:

    mutable ListenerList m_Listeners;
    mutable Poco::FastMutex m_Mutex;

};

// message with 4 parameters and return type
template <typename T, typename U, typename V, typename W, typename A = void>
class Message4
{
  public:

    typedef Message4 Self;
    typedef MessageAbstractDelegate4<T,U,V,W,A> AbstractDelegate;
    typedef std::vector<AbstractDelegate* > ListenerList;

    ~Message4() {
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end(); ++iter )
      {
        delete *iter;
      }
    }

    void AddListener( const AbstractDelegate& delegate ) const
    {
      AbstractDelegate* msgCmd = delegate.Clone();

      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(msgCmd)) {
          delete msgCmd;
          return;
        }
      }
      m_Listeners.push_back(msgCmd);
    }

    void operator += ( const AbstractDelegate& delegate ) const
    {
      this->AddListener(delegate);
    }

    void RemoveListener( const AbstractDelegate& delegate ) const
    {
      Poco::FastMutex::ScopedLock lock(m_Mutex);
      for (typename ListenerList::iterator iter = m_Listeners.begin();
           iter != m_Listeners.end();
           ++iter )
      {
        if ((*iter)->operator==(&delegate))
        {
          delete *iter;
          m_Listeners.erase( iter );
          return;
        }
      }
    }

    void operator -= ( const AbstractDelegate& delegate) const
    {
      this->RemoveListener(delegate);
    }

    void Send(T t, U u, V v, W w)
    {
      ListenerList listeners;

      {
        Poco::FastMutex::ScopedLock lock(m_Mutex);
        listeners.assign(m_Listeners.begin(), m_Listeners.end());
      }

      for ( typename ListenerList::iterator iter = listeners.begin();
            iter != listeners.end();
            ++iter )
      {
        // notify each listener
        (*iter)->Execute(t,u,v,w);
      }
    }

    void operator()(T t, U u, V v, W w)
    {
      this->Send(t, u, v, w);
    }

    const ListenerList& GetListeners()
    {
      return this->m_Listeners;
    }

  protected:

    mutable ListenerList m_Listeners;
    mutable Poco::FastMutex m_Mutex;

};

} // namespace cherry

#endif /*CHERRYMESSAGE_H_*/
